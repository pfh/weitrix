

# Least squares and then apply linear combinations, 
# calculate vcov.
# K %*% weighted_pseudoinverse(A) %*% b
fit_and_contrast <- function(K, X, w, y) {
    present <- w > 0
    w <- w[present]
    y <- y[present]
    X <- X[present,,drop=FALSE]

    # Insufficient data?
    if (length(y) < ncol(X)) 
        return(NULL)

    sw <- sqrt(w)
    swy <- sw*y
    swX <- sw*X
    decomp <- svd(swX)
    
    # Fail if singular
    if (any(decomp$d == 0)) 
        return(NULL) 

    solver <- decomp$v %*% (t(decomp$u)/decomp$d)
    
    beta_hat <- as.vector(solver %*% swy)
    residuals <- swy - swX %*% beta_hat      #(weighted)
    estimates <- as.vector(K %*% beta_hat)

    Ksolver <- K %*% solver
    unscaled_covar <- Ksolver %*% t(Ksolver)

    # How big is the linear combination?
    eig <- eigen(unscaled_covar)
    unscaled_ss <- sum( (t(eig$vectors) %*% estimates)^2 / eig$values )
    # unscaled_ss <- t(estimates) %*% solve(unscaled_cover) %*% estimates

    # A unit step could cause up to this much inflation in SS
    #biggest_unit_step <- 1/min(eig$values)

    list(
        estimates = estimates,
        unscaled_vars = diag(unscaled_covar),
        #unscaled_covar = unscaled_covar,
        eigenvalues = eig$values,
        unscaled_ss = unscaled_ss,
        sum_weight = sum(w),
        #biggest_unit_step = biggest_unit_step,
        #kappa = max(eig$value)/min(eig$values),
        n_present = length(y),
        df_residuals = length(y)-ncol(X),
        ss_residuals = sum(residuals^2),
        weighted_mean = weighted.mean(y,w)
    )
}

fit_and_contrast_inner <- function(args) with(args, {
    n <- nrow(Y)
    result <- rep(list(NULL), n)
    for(i in seq_len(n))
        result[[i]] <- fit_and_contrast(K,X,W[i,],Y[i,])
    result
})


fit_and_contrast_all <- function(K, X, W, Y) {
    this_function_bp_up()
    BPPARAM <- bpparam()

    parts <- partitions(nrow(Y), ncol(Y)*2, cpu_heavy=TRUE, BPPARAM=BPPARAM)
    feed <- map(parts, function(part) {
        list(
            K=K, X=X,
            W=W[part,,drop=FALSE],
            Y=Y[part,,drop=FALSE])
    })

    result <- bplapply(feed, fit_and_contrast_inner, BPPARAM=BPPARAM)
    do.call(c, result)
}


mvtreat_inner <- function(q, ss_observed, ss_hypothesized, df1, df2) {
    precision <- qchisq(q, df=df2)/df2
    1 - pchisq(ss_observed*precision, df=df1, ncp=ss_hypothesized*precision)
}

# Multivariate extension of TREAT
#
# This is a weird fusion of Bayesian and Frequentist
# - Bayesian in that we have a prior belief about the distribution of residual variance
#   which should be updated by observing the actual residuals before calling
#   this function. (Typically this will be from empirical Bayes a la limma.)
# - Frequentist in that we have no distributional beliefs about sum of squares of
#   the quantities we're estimating.
#
# We believe the variance to be drawn from an inverse gamma distribution, 
# as if produced by var*df2/rchisq(1,df2)
#
# We hypothesize the sum of squares to be less than or equal to ss_hypothesized*var
# We observe a sum of squares ss_observed*var
#
# (Note: divide out "var" from everything before calling.)
#
pmvtreat <- function(ss_observed, ss_hypothesized, df1, df2) {
    if (ss_hypothesized == 0)
        return(pf(ss_observed/df1, df1, df2, lower.tail=FALSE))
    
    if (ss_observed <= ss_hypothesized)
        return(1)

    if (df2 == Inf)
        return(1 - pchisq(ss_observed, df=df1, ncp=ss_hypothesized))
    
    result <- integrate(mvtreat_inner, 0, 1, 
        ss_observed, ss_hypothesized, df1, df2, 
        abs.tol=0, stop.on.error=FALSE)

    #if (result$message != "OK")
    #    warning(paste0("While calulating p-value, ", result$message))

    result$value
}

#' Top confident effects based on one or more contrasts of a linear model for each row
#'
#' This function provides topconfects-style testing of a linear model contrast, as well as a multi-contrast extension of this method for F-tests with effect sizes.
#'
#' Based on the \code{effect} argument, the estimated effect may be:
#'
#' \itemize{
#' \item{\code{"auto"}} Choose \code{"contrast"} or \code{"sd"} as appropriate. 
#' \item{\code{"contrast"}} The estimated contrast. This should produce results identical to a limma-topconfects analysis.
#' \item{\code{"sd"}} Standard deviation explained (i.e. square root of the variance explained) by the part of the model captured by the contrasts provided. 
#' \item{\code{"cohen_f"}} Cohen's f, i.e. the signal to noise ratio. Ranking is similar to traditional ranking of results by p-value. 
#' }
#'
#' Based on the \code{dispersion_est} argument, the estimated residual dispersion is estimated as:
#'
#' \itemize{
#' \item{\code{"none"}} Weitrix is assumed to be fully calibrated already. Dispersion is assumed to be 1. If the assumption is correct, this is most powerful, as there is no uncertainty to the dispersion.
#' \item{\code{"row"}} Dispersion is estimated based on the residuals for each row. With a limited number of columns, this estimate is uncertain (low residual degrees of freedom), so may lack power.
#' \item{\code{"ebayes_limma"}} Default, recommended. Perform Empricial Bayes squeezing of dispersions, using \code{limma::squeezeVar}. This also reduces the uncertainty about the dispersion (mainfesting as extra "prior" degrees of freedom), increasing the power of the test.
#' }
#'
#' In results from this function, whenever we talk about the mean, standard deviation explained, or typical observation error, this should be understood to be weighted. Standard deviation explained is in the same units as the observations, but its estimation is weighted by the weights, so in a row with some high weight observations and other low weight observations, estimated standard deviation explained will mostly be driven by the high weight observations.
#'
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#' with \code{as_weitrix}.
#' @param design 
#' A formula in terms of \code{colData(weitrix} or a design matrix, 
#' which will be fitted to the weitrix on each row. 
#' Can also be a pre-existing Components object, 
#' in which case the existing fits (\code{design$row}) are used.
#' @param coef
#' Give either coef or contrasts but not both.
#' If coef is given, it is converted into a set of contrasts that simply
#'   test each given coefficient. 
#' Coefficients can be specified by number of name.
#' @param contrasts
#' Give either coef or contrasts but not both.
#' One or more contrasts of interest, i.e. specifications of linear combination of coefficients.
#' Each contrast should be placed in a columns.
#' The number of rows should match the number of coefficients.
#' @param effect
#' Effect to estimate and provide confidence bounds on. See description.
#' @param dispersion_est
#' Method of estimating per-row dispersion. See description.
#' @param fdr 
#' False Discovery Rate to control for.
#' @param step 
#' Granularity of effect sizes to test.
#' @param full 
#' If TRUE, output some further columns related to the calculations. 
#'
#' @return
#' A topconfects result. The \code{$table} data frame contains columns:
#' \itemize{
#' \item{effect}{ Estimated effect (as requested using the \code{effect} parameter). }
#' \item{confect}{ An inner confidence bound on effect. }
#' \item{fdr_zero}{ FDR-adjusted p-value for the null hypothesis that effect is zero. }
#' \item{row_mean}{ Weighted row mean. }
#' \item{typical_obs_err}{ Typical residual standard deviation (square root of variance) associated with observations in this row. Note that each observation has its own associated variance, based on its weight and the row dispersion estimate used. This column is calculated from the weighted average variance of observations. }
#' } 
#'
#' @examples
#'
#' # Simplest possible test
#' # Which rows have an average different from zero?
#' weitrix_confects(simwei, ~1, coef="(Intercept)")
#'
#' # See vignettes for more substantial examples
#'
#' @export
weitrix_confects <- function(
        weitrix, design, coef=NULL, contrasts=NULL, 
        effect=c("auto","contrast","sd","cohen_f"),
        dispersion_est=c("ebayes_limma","row","none"),
        fdr=0.05, step=NULL, full=FALSE) {

    effect_wanted <- match.arg(effect)
    dispersion_est <- match.arg(dispersion_est)
    assert_that(!is.null(coef) || !is.null(contrasts))
    assert_that(is.null(coef) || is.null(contrasts))

    weitrix <- as_weitrix(weitrix)

    if (is(design, "formula"))
        design <- model.matrix(design, data=colData(weitrix))

    # Convert coefficients into contrast matrix if given
    if (!is.null(coef)) {
        mat <- diag(ncol(design))
        colnames(mat) <- colnames(design)
        contrasts <- mat[, coef, drop=FALSE]
    }

    design <- as.matrix(design)
    contrasts <- as.matrix(contrasts)
    n <- nrow(weitrix)

    assert_that(ncol(weitrix) == nrow(design))
    assert_that(nrow(contrasts) == ncol(design))

    if (effect_wanted == "auto") {
        if (ncol(contrasts) == 1)
            effect_wanted <- "contrast"
        else
            effect_wanted <- "sd"
    }

    if (is.null(step)) {
        if (effect_wanted == "contrast")
            step <- 0.001
        else
            step <- 0.01
    }

    assert_that(effect_wanted != "contrast" || ncol(contrasts) == 1, 
        msg='Only one contrast allowed for "contrast" effect.') 


    # ==== Fit and contrast each row ====

    fits <- fit_and_contrast_all( 
        t(contrasts), design, weitrix_weights(weitrix), weitrix_x(weitrix)) 


    # ==== Extract vectors from list of fits ====
    
    good <- !map_lgl(fits, is.null)
    fits_good <- fits[good]
    get <- function(what) {
        x <- rep(NA_real_, n)
        x[good] <- map_dbl(fits_good, what)
        x
    }

    ss1 <- get("unscaled_ss")
    df1 <- ncol(contrasts)
    sum_weight <- get("sum_weight")
    n_present <- get("n_present")
    ss_residuals <- get("ss_residuals")
    df_residuals <- get("df_residuals")


    # ==== Define denominator of F ratio as requested ====

    if (dispersion_est == "none") {
        sigma2 <- rep(1, n)
        df2 <- rep(Inf, n)
    } else if (dispersion_est == "row") {
        sigma2 <- ss_residuals / df_residuals
        df2 <- df_residuals
        sigma2[ df2==0 ] <- NA
    } else if (dispersion_est == "ebayes_limma") {
        # Which rows can we actually estimate the residual variance for?
        very_good <- good & df_residuals > 0
        assert_that(sum(very_good) > 0, msg="No remaining degrees of freedom.")
        squeeze <- squeezeVar(
            ss_residuals[very_good]/df_residuals[very_good], 
            df_residuals[very_good])
        df2 <- df_residuals + squeeze$df.prior
        if (squeeze$df.prior == Inf)
            sigma2 <- rep(squeeze$var.prior, n)
        else
            sigma2 <- (ss_residuals+squeeze$var.prior*squeeze$df.prior)/df2
    }


    # ==== Calculate F ratio ====

    F <- (ss1/df1)/sigma2


    if (effect_wanted == "contrast") {
       # ==== Topconfects for single contrast =====

       effect_desc = "contrast"
       effect <- get("estimates")
       se <- sqrt( get("unscaled_vars")*sigma2 )

       result <- normal_confects(
           effect, se, df=df2, fdr=fdr, step=step, full=TRUE)
    } else if (effect_wanted == "cohen_f") {
        # ==== Topconfects for one or more contrasts ====

        effect_desc <- "Cohen's f"
        effect <- sqrt( (F*df1)/n_present )

        pfunc <- function(indices, mag) {
            1 - pf(
                F[indices],
                df1=df1,
                df2=df2[indices],
                ncp=mag^2 * n_present[indices])
        }

        result <- nest_confects(
            n, pfunc, fdr=fdr, step=step, full=TRUE)
    } else {
        # ==== Mutlivariate version of TREAT ====

        effect_desc <- "standard deviation explained"
        effect <- sqrt( ss1/sum_weight )

        pfunc <- function(indices, mag) {
            map_dbl(indices, ~pmvtreat(
                ss_observed=ss1[.]/sigma2[.],
                ss_hypothesized=mag^2*sum_weight[.]/sigma2[.],
                df1=df1,
                df2=df2[.]))
        }

        result <- nest_confects(
            n, pfunc, fdr=fdr, step=step, full=TRUE)
    }


    # ==== Provide further information in result ====
    
    fdr_zero <- result$table$fdr_zero
    result$table$fdr_zero <- NULL

    result$table$effect <- effect[result$table$index]

    if (full || effect_wanted != "contrast") {
        for(i in seq_len(ncol(contrasts))) {
            name <- colnames(contrasts)[i]
            if (is.null(name)) name <- paste0("contrast",i)
            result$table[[name]] <- get(~.$estimates[i])[result$table$index]
        }
    }

    if (full) {
        for(i in seq_len(ncol(contrasts))) {
            name <- colnames(contrasts)[i]
            if (is.null(name)) name <- paste0("contrast",i)
            result$table[[paste0(name,"_se")]] <- 
                sqrt(get(~.$unscaled_vars[i])*sigma2)[result$table$index]
        }
    }

    result$table$fdr_zero <- fdr_zero
    result$table$row_mean <- get("weighted_mean")[result$table$index]
    result$table$typical_obs_err <- sqrt(sigma2*n_present/sum_weight)[result$table$index]
    if (full) {
        result$table$F <- F[result$table$index]
        result$table$n_present <- n_present[result$table$index]
        result$table$dispersion_seen <- (ss_residuals/df_residuals)[result$table$index]
        result$table$df_seen <- df_residuals[result$table$index]
        result$table$dispersion_used <- sigma2[result$table$index]
        result$table$df_used <- df2[result$table$index]
    }

    result$table$name <- rownames(weitrix)[result$table$index]

    for(name in colnames(rowData(weitrix)))
        if (!name %in% colnames(result$table))
            result$table[[name]] <- rowData(weitrix)[result$table$index,name]

    result$effect_desc <- effect_desc

    # TODO: report these details when printed
    result$dispersion_est <- dispersion_est
    result$df_top <- df1
    if (dispersion_est == "ebayes_limma") {
        result$df_prior <- squeeze$df.prior
        result$var_prior <- squeeze$var.prior
    }

    result$design <- design
    result$contrasts <- contrasts

    result
}



