
# Weighted variance using a matrix decomposition as model
calc_row_dispersion_inner <- function(args) with(args, {
    x <- as.matrix(x)
    w <- as.matrix(w)
    n <- nrow(x)
    p <- ncol(row)

    result <- rep(NA_real_,n)
    for(i in seq_len(n)) {
        wi <- w[i,]
        present <- wi > 0
        wi_present <- wi[present]
        df <- length(wi_present)-p
        if (df > 0) {
            errors <- as.vector(x[i,present]) - 
                as.vector(col[present,,drop=FALSE] %*% row[i,]) 
            result[i] <- sum(errors*errors*wi_present) / df
        }
    }
    result
})

calc_row_dispersion <- function(x, w, row, col) {
    BPPARAM <- bpparam()

    parts <- partitions(nrow(x), ncol(x)*2, BPPARAM=BPPARAM)
    feed <- map(parts, function(part) {
        list(
            x=x[part,,drop=FALSE],
            w=w[part,,drop=FALSE],
            row=row[part,,drop=FALSE],
            col=col)
    })
    result <- bplapply(feed, calc_row_dispersion_inner, BPPARAM=BPPARAM)
    unlist(result)
}


as_components <- function(design, weitrix) {
    if (is(design, "Components"))
        result <- design
    else
        result <- weitrix_components(weitrix, design=design, p=0, verbose=FALSE)
    
    assert_that(nrow(result$row) == nrow(weitrix))
    assert_that(nrow(result$col) == ncol(weitrix))
    
    result
}


#' Calculate row dispersions
#'
#' Calculate the dispersion of each row. For each observation, 
#' this value divided by the weight gives the observation's variance.
#'
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#' with \code{as_weitrix}.
#' @param design 
#' A formula in terms of \code{colData(weitrix} or a design matrix, 
#' which will be fitted to the weitrix on each row. 
#' Can also be a pre-existing Components object, 
#' in which case the existing fits (\code{design$row}) are used.
#'
#' @return
#' A numeric vector.
#'
#' @examples
#' # Using a model just containing an intercept
#' weitrix_dispersions(simwei, ~1)
#'
#' # Allowing for one component of variation, the dispersions are lower
#' comp <- weitrix_components(simwei, p=1, verbose=FALSE)
#' weitrix_dispersions(simwei, comp)
#'
#' @export
weitrix_dispersions <- function(weitrix, design=~1) {
    weitrix <- as_weitrix(weitrix)
    comp <- as_components(design, weitrix)

    calc_row_dispersion(
        weitrix_x(weitrix),
        weitrix_weights(weitrix),
        comp$row,
        comp$col)        
}


#' Adjust weights row-wise based on given row dispersions
#'
#' Based on estimated row dispersions, adjust weights in each row.
#'
#' For large numbers of samples this can be based directly 
#' on weitrix_dispersions. 
#' For small numbers of samples, when using limma, 
#' it should be based on a trend-line fitted to known co-variates of 
#' the dispersions. 
#' This can be done using \code{weitrix_calibrate_trend}.
#'
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#' with \code{as_weitrix}.
#' @param dispersions 
#' A dispersion for each row.
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#' 
#' @examples
#' # Adjust weights so dispersion for each row is exactly 1. This is dubious 
#' # for a small dataset, but would be fine for a dataset with many columns.
#' comp <- weitrix_components(simwei, p=1, verbose=FALSE)
#' disp <- weitrix_dispersions(simwei, comp)
#' cal <- weitrix_calibrate(simwei, disp)
#' weitrix_dispersions(cal, comp)
#'
#' @export
weitrix_calibrate <- function(weitrix, dispersions) {
    weitrix <- as_weitrix(weitrix)
    assert_that(nrow(weitrix) == length(dispersions))

    # Zero out weights if dispersion not available
    mul <- 1/dispersions
    mul[is.na(mul)] <- 0

    weitrix_weights(weitrix) <- sweep(
        weitrix_weights(weitrix), 1, mul, "*")
    
    metadata(weitrix)$weitrix$calibrated <- TRUE
    
    weitrix
}




#' Adjust weights row-wise by fitting a trend to estimated dispersions
#'
#' Dispersions are estimated using \code{weitrix_dispersions}. 
#' A trend line is then fitted to the dispersions using a gamma GLM
#'     with log link function.
#' Weitrix weights are calibrated based on this trend line. 
#'
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#' with \code{as_weitrix}.
#' @param design 
#' A formula in terms of \code{colData(weitrix} or a design matrix, 
#' which will be fitted to the weitrix on each row. 
#' Can also be a pre-existing Components object, 
#' in which case the existing fits (\code{design$row}) are used.
#' @param trend_formula 
#' A formula specification for predicting log dispersion from 
#' columns of rowData(weitrix). 
#' If absent, metadata(weitrix)$weitrix$calibrate_trend_formula is used.
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#' 
#' Several columns are added to the \code{rowData}:
#' \itemize{
#'   \item{deg_free}{ Degrees of freedom for dispersion calculation.}
#'   \item{dispersion_before}{ Dispersion before calibration.}
#'   \item{dispersion_trend}{ Fitted dispersion trend.}
#'   \item{dispersion_after}{ Dispersion for these new weights.}
#' }
#'
#' @examples
#' rowData(simwei)$total_weight <- rowSums(weitrix_weights(simwei))
#'
#' # To estimate dispersions, use a simple model containing only an intercept
#' # term. Model log dispersion as a straight line relationship with log total 
#' # weight and adjust weights to remove any trend. 
#' cal <- weitrix_calibrate_trend(simwei,~1,trend_formula=~log(total_weight))
#'
#' # This dataset has few rows, so calibration like this is dubious.
#' # Predictors in the fitted model are not significant.
#' summary( metadata(cal)$weitrix$trend_fit )
#'
#' # Information about the calibration is added to rowData
#' rowData(cal)
#'
#'
#' # A Components object may also be used as the design argument.
#' comp <- weitrix_components(simwei, p=1, verbose=FALSE)
#' cal2 <- weitrix_calibrate_trend(simwei,comp,trend_formula=~log(total_weight))
#'
#' rowData(cal2)
#'
#' @export
weitrix_calibrate_trend <- function(weitrix, design=~1, trend_formula=NULL) {
    weitrix <- as_weitrix(weitrix)
    comp <- as_components(design, weitrix)

    if (is.null(trend_formula))
        trend_formula <- metadata(weitrix)$weitrix$calibrate_trend_formula
    if (is.null(trend_formula))
        trend_formula <- metadata(weitrix)$weitrix$trend_formula
    assert_that(!is.null(trend_formula))
    
    trend_formula <- as.formula(trend_formula)

    rowData(weitrix)$deg_free <- 
        rowSums(weitrix_weights(weitrix) > 0) - ncol(comp$col)
    rowData(weitrix)$dispersion_before <- weitrix_dispersions(weitrix, comp)
    
    data <- rowData(weitrix)
    
    # Least squares method on log dispersion_before
    #
    #trend_formula <- update(trend_formula, log(dispersion_before)~.)
    #
    #tiny <- mean(data$dispersion_before,na.rm=TRUE)*1e-9
    #good <- !is.na(data$dispersion_before) & data$dispersion_before > tiny
    #data <- data[good,]
    #
    #fit <- eval(substitute(
    #    lm(trend_formula, data=data),
    #    list(trend_formula=trend_formula)
    #))
    #
    #pred <- exp(
    #    predict(fit, newdata=rowData(weitrix)) +
    #    sigma(fit)^2/2)
    

    # Gamma GLM method
    # - Gamma() doesn't like zeros, but it's actually fine, 
    #   so use quasi(), fit will be identical.
    # - Choice of starting predictions for iteration matters,
    #   very small dispersions cause numerical errors,
    #   so start with a conservative set of predictions. 

    trend_formula <- update(trend_formula, dispersion_before~.)
    mustart <- mean(data$dispersion_before, na.rm=TRUE)
    n <- nrow(data)

    fit <- eval(substitute(
        glm2(
            trend_formula, 
            data=data,
            weights=deg_free,
            family=quasi(link="log", variance="mu^2"),
            mustart=rep(mustart, n),
            control=glm.control(maxit=100)
        ),
        list(trend_formula=trend_formula, mustart=mustart, n=n)
    ))

    pred <- predict(fit, newdata=rowData(weitrix), type="response")

    rowData(weitrix)$dispersion_trend <- pred
    rowData(weitrix)$dispersion_after <- 
        rowData(weitrix)$dispersion_before / pred
    
    metadata(weitrix)$weitrix$trend_fit <- fit
    
    weitrix_calibrate(weitrix, pred)
}


# Hack!
# Extract anything that might be referred to in a formula
all_names <- function(obj) {
    if (inherits(obj, "name"))
        return(as.character(obj))
    
    if (!inherits(obj, c("formula","call")))
        return(c())
    
    do.call(c, lapply(obj, all_names))
}


#' Adjust weights element-wise by fitting a trend to squared residuals
#'
#' This is a very flexible method of calibrating weights.
#' It should be especially useful if your existing weights account for 
#'     technical variation, 
#'     but there is also biological variation.
#' In this case large weights will tend to be overly optimistic, 
#'     and a non-linear transformation of weights is needed.
#'
#' Residuals are found relative to a fitted model.
#' A trend model is then fitted to the squared residuals using a gamma GLM
#'     with log link function.
#' Weitrix weights are set based on the inverse of the fitted trend. 
#'
#' Residuals from a fitted model are generally smaller than residuals 
#'     from the true model.
#' A simple adjustment to the weights is made to account for this.
#' Weights are reduced by a factor of (n-ncol(design)*nrow(weitrix))/n 
#'     where n is the number of non-missing values in the weitrix.
#'
#' \code{trend_formula} may reference any row or column variables,
#'     or \code{mu} for the predicted value,
#'     or \code{weight} for the existing weights,
#'     or special factors \code{row} and \code{col}.
#' Keep in mind also that a log link function is used.
#' 
#' Unlike in \code{weitrix_calibrate_trend},
#'     existing weights must be explicitly included in the formula
#'     if they are to be retained (see examples).
#' 
#' This function is currently not memory efficient,
#'     it should be fine for bulk experiments but may struggle for single cell.
#' To reduce memory usage somewhat, 
#'     when constructing the data frame on which to fit the glm, 
#'     only columns referenced in \code{trend_formula} are included.
#'
#' Example formulas:
#'
#' \code{trend_formula=~1+offset(-log(weight))} 
#' Apply a global scaling, otherwise keeping weights the same.
#'
#' \code{trend_formula=~log(weight)}
#' Moderate weights by raising them to some power 
#' and applying some overall scaling factor. 
#' This will allow for biological variation.
#'
#' \code{trend_formula=~poly(log(weight),2))}
#' Apply a more complex quadratic curve-based moderation of weights.
#'
#' \code{trend_formula=~col+offset(-log(weight))}
#' Calibrate each sample's weights by a scaling factor.
#' Note that due to the simplistic adjustment for using a fitted model rather
#'     than the true model, this may give misleading results when 
#'     the design is unbalanced and there are few samples, 
#'     i.e. when there are some samples with much higher leverage than others.
#'
#' \code{trend_formula=~col*poly(log(weight),2)}
#' Quadratic curve moderation of weights, applied to each sample individually.
#'
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#' with \code{as_weitrix}.
#' @param design 
#' A formula in terms of \code{colData(weitrix} or a design matrix, 
#' which will be fitted to the weitrix on each row. 
#' Can also be a pre-existing Components object, 
#' in which case the existing fits (\code{design$row}) are used.
#' @param trend_formula 
#' A formula specification for predicting squared residuals. See below.
#' If absent, metadata(weitrix)$weitrix$calibrate_all_formula is used.
#'
#' @param mu_min
#' When fitting the GLM, omit observations 
#'     where the estimated mu is less than this value.
#' When calculating weights from the fitted GLM, 
#'     clip mu to be at least this value.
#' @param mu_max
#' When fitting the GLM, omit observations 
#'     where the estimated mu is greater than this value.
#' When calculating weights from the fitted GLM, 
#'     clip mu to be at most this value.
#'
#' @param keep_fit
#' Keep glm fit and the data used to create it. This can be large!
#' If TRUE, these will be stored in \code{metadata(weitrix)$weitrix$all_fit}
#' and \code{metadata(weitrix)$weitrix$all_data}.
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#'
#' \code{metadata(weitrix)$weitrix} will contain the fitted trend model,
#'    and if requested the data frame used to fit the model.
#'
#' @examples
#'
#' simcal <- weitrix_calibrate_all(simwei, ~1, ~log(weight), keep_fit=TRUE)
#'
#' metadata(simcal)$weitrix$all_fit
#'
#' @export
weitrix_calibrate_all <- function(
        weitrix, design=~1, trend_formula=NULL, 
        mu_min=NA, mu_max=NA,
        keep_fit=FALSE) {

    weitrix <- as_weitrix(weitrix)
    comp <- as_components(design, weitrix)

    if (is.null(trend_formula))
        trend_formula <- metadata(weitrix)$weitrix$calibrate_all_formula
    assert_that(!is.null(trend_formula))
    
    trend_formula <- as.formula(trend_formula)

    needed <- all_names(trend_formula)
    needed_rowdata <- intersect(colnames(rowData(weitrix)), needed)
    needed_coldata <- intersect(colnames(colData(weitrix)), needed)

    data <- cbind(
        as.data.frame(rowData(weitrix)[
            rep(seq_len(nrow(weitrix)), ncol(weitrix)),
            needed_rowdata,
            drop=FALSE]),
        as.data.frame(colData(weitrix)[
            rep(seq_len(ncol(weitrix)), each=nrow(weitrix)),
            needed_coldata,
            drop=FALSE])
    )

    if ("row" %in% needed)
        data$row <- rep(
            factor(rownames(weitrix), rownames(weitrix)), ncol(weitrix))

    if ("col" %in% needed)
        data$col <- rep(
            factor(colnames(weitrix), colnames(weitrix)), each=nrow(weitrix))
    
    data$weight <- as.vector(weitrix_weights(weitrix))
    data$mu <- as.vector(comp$row %*% t(comp$col))
    data$.y <- (as.vector(weitrix_x(weitrix)) - data$mu)^2

    # mu may be clipped to a range of values
    # - data outside the range is treated as missing
    # - when assigning weights, mu is clipped to this range
    if (!is.na(mu_min)) {
        clip <- data$mu < mu_min
        data$.y[clip] <- NA
        data$mu[clip] <- mu_min
    }

    if (!is.na(mu_max)) {
        clip <- data$mu > mu_max
        data$.y[clip] <- NA
        data$mu[clip] <- mu_max
    }

    # Want to be able to use poly(log(weight),2)
    # so zero weights need to be dropped.
    good <- data$weight > 0 
    good_data <- data[good,,drop=FALSE]

    mustart <- mean(data$.y, na.rm=TRUE)
    n <- nrow(good_data)

    trend_formula <- update(trend_formula, .y ~ .)
    
    fit <- eval(substitute(
        glm2(
            trend_formula,
            data=good_data,
            family=quasi(link="log", variance="mu^2"),
            mustart=rep(mustart, n),
            model=FALSE, x=FALSE, y=FALSE,
            control=glm.control(maxit=100)
        ),
        list(trend_formula=trend_formula, mustart=mustart, n=n)
    ))

    n_good <- nrow(good_data)
    overfitting_adjustment <- (n_good-ncol(comp$col)*nrow(weitrix)) / n_good

    pred <- predict(fit, newdata=good_data, type="response")
    pred[is.na(pred)] <- Inf
    new_weights <- matrix(0, nrow=nrow(weitrix), ncol=ncol(weitrix))
    new_weights[good] <- overfitting_adjustment/pred
    rownames(new_weights) <- rownames(weitrix)
    colnames(new_weights) <- colnames(weitrix)

    weitrix_weights(weitrix) <- new_weights

    metadata(weitrix)$weitrix$all_coef <- coef(fit)

    if (keep_fit) {
        metadata(weitrix)$weitrix$all_fit <- fit
        data$new_weight <- as.vector(weitrix_weights(weitrix))
        metadata(weitrix)$weitrix$all_data <- data
    }

    weitrix
}


#' Weight calibration plots, optionally versus a covariate
#'
#' Various plots based on weighted squared residuals 
#'     of each element in the weitrix.
#' \code{weight*residual^2} is the Pearson residual for a gamma GLM plus one,
#'     as used by \code{weitrix_calibrate_all}.
#'
#' This function is not memory efficient.
#' It is suitable for typical bulk data, but generally not not for single-cell.
#'
#' Defaults to a boxplot of sqrt(weight) weighted residuals.
#' Blue guide bars are shown for the expected quartiles,
#'     these will ideally line up with the boxplot.
#'
#' If \code{cat} is given, it will be used to break the elements down
#'     into categories.
#'
#' If \code{covar} is given, sqrt(weight) weighted residuals are plotted versus
#'     the covariate, with red trend lines for the mean and for the mean +/- one standard deviation.
#' If the weitrix is calibrated, the trend lines should be horizontal lines
#'     with y intercept close to -1, 0 and 1.
#' Blue guide lines are shown for this ideal outcome.
#'
#' Any of the variables available with \code{weitrix_calibrate_all}
#'     can be used for \code{covar} or \code{cat}.
#'
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#' with \code{as_weitrix}.
#' @param design 
#' A formula in terms of \code{colData(weitrix} or a design matrix, 
#' which will be fitted to the weitrix on each row. 
#' Can also be a Components object.
#' @param covar
#' Optional. A covariate. Specify as you would with \code{ggplot2::aes}.
#' Can be a matrix of the same size as \code{weitrix}.
#' @param cat
#' Optional. A categorical variable to break down the data by.
#' Specify as you would with \code{ggplot2::aes}.
#' @param funnel
#' Flag. Produce a funnel plot? 
#' Note: \code{covar} can not be used for funnel plots.
#' @param guides
#' Show blue guide lines.
#'
#' @return
#' A ggplot2 plot.
#'
#' @examples
#' weitrix_calplot(simwei, ~1)
#' weitrix_calplot(simwei, ~1, covar=mu)
#' weitrix_calplot(simwei, ~1, cat=col)
#'
#' # weitrix_calplot should generally be used after calibration
#' cal <- weitrix_calibrate_all(simwei, ~1, ~col+log(weight))
#' weitrix_calplot(cal, ~1, cat=col)
#'
#' # You can use a matrix of the same size as the weitrix as a covariate.
#' # It will often be useful to assess vs the original weighting.
#' weitrix_calplot(cal, ~1, covar=weitrix_weights(simwei))
#'
#' @export
weitrix_calplot <- function(
        weitrix, design=~1, covar, cat, funnel=FALSE, guides=TRUE) {

    covar_var <- enquo(covar)
    cat_var <- enquo(cat)
    have_covar <- !identical(covar_var, quo())
    have_cat <- !identical(cat_var, quo())

    if (funnel) {
        assert_that(!have_covar, msg="Can't use covar with funnel.")
        weight <- NULL #Prevent R CMD check note. Is there a better way?
        covar_var <- quo(1/sqrt(weight))
        have_covar <- TRUE
    }

    weitrix <- as_weitrix(weitrix)
    comp <- as_components(design, weitrix)

    needed <- c(all_names(covar_var), all_names(cat_var))
    needed_rowdata <- intersect(colnames(rowData(weitrix)), needed)
    needed_coldata <- intersect(colnames(colData(weitrix)), needed)

    data <- cbind(
        as.data.frame(rowData(weitrix)[
            rep(seq_len(nrow(weitrix)), ncol(weitrix)),
            needed_rowdata,
            drop=FALSE]),
        as.data.frame(colData(weitrix)[
            rep(seq_len(ncol(weitrix)), each=nrow(weitrix)),
            needed_coldata,
            drop=FALSE])
    )

    if ("row" %in% needed)
        data$row <- 
            rep(factor(rownames(weitrix), rownames(weitrix)), ncol(weitrix))
    
    if ("col" %in% needed)
        data$col <- 
            rep(factor(colnames(weitrix), colnames(weitrix)), each=nrow(weitrix))

    data$weight <- as.vector(weitrix_weights(weitrix))
    data$mu <- as.vector(comp$row %*% t(comp$col))
    data$residual <- ifelse(
        data$weight > 0,
        as.vector(weitrix_x(weitrix)) - data$mu,
        rep(NA, nrow(data)))
    data$weighted_residual <- sqrt(data$weight) * data$residual

    n_good <- sum(!is.na(data$residual))
    overfitting_adjustment <- (n_good-ncol(comp$col)*nrow(weitrix)) / n_good
    guide_pos <- sqrt(overfitting_adjustment)
    
    if (!have_cat) {
        data$weitrix <- ""
        cat_var <- quo(weitrix)
    }

    if (have_covar) {
        y <- if (funnel) data$residual else data$weighted_residual
        x <- eval_tidy(covar_var, data=data)
        cat <- droplevels(factor(eval_tidy(cat_var, data=data)))
        present <- !is.na(y)
        y <- y[present]
        x <- x[present]
        cat <- cat[present]
        plot_data <- data.frame(x=x,y=y,cat=cat)

        # Calculate a trend line, needs to be square-unbiassed
        # Could be made more efficient
        bin_cat <- c()
        bin_x <- c()
        bin_mean <- c()
        bin_sd <- c()
        for(level in levels(cat)) {
            this_x <- x[ cat == level ]
            this_y <- y[ cat == level ]
            n <- length(this_x)
            n_bins_wanted <- ceiling(n^(1/3))
            bins <- droplevels(cut(rank(this_x), 
                seq(0.5,n+0.5,length.out=n_bins_wanted+1)))
            bin_cat <- c(bin_cat, rep(level, length(levels(bins))))
            bin_x <- c(bin_x, tapply(this_x, list(bins), mean))
            bin_mean <- c(bin_mean, tapply(this_y, list(bins), mean))
            bin_sd <- c(bin_sd, tapply(this_y, list(bins), sd))
        }
        bin_cat <- factor(bin_cat, levels=levels(cat))

        line_data <- data.frame(x=bin_x,mean=bin_mean,sd=bin_sd,cat=bin_cat)

        ggplot(plot_data, aes(y=.data$y, x=.data$x)) + 
            {if (have_cat) facet_wrap(vars(.data$cat))} +
            #geom_point(stroke=0, size=0.5, alpha=0.5, na.rm=TRUE) + 
            geom_bin2d(bins=200) + scale_fill_viridis_c() +
            {if (guides && !funnel) geom_hline(yintercept=c(-guide_pos,0,guide_pos), color="blue")} + 
            {if (guides && funnel) geom_abline(slope=c(guide_pos,0,-guide_pos),intercept=c(0,0,0),color="blue")} +
            geom_line(aes(y=.data$mean), data=line_data, size=1, color="red")+
            geom_line(aes(y=.data$mean+.data$sd), data=line_data, size=1, color="red")+
            geom_line(aes(y=.data$mean-.data$sd), data=line_data, size=1, color="red")+
            labs(
                y=if (funnel) "residual" else "sqrt(weight) * residual", 
                x=as_label(covar_var),
                fill="Count")+
            theme(legend.position="bottom", legend.justification=c(1,1))
    } else {
        ggplot(data, aes(
                y=.data$weighted_residual, 
                x=!!cat_var)) + 
            {if (guides) geom_hline(
                yintercept=guide_pos*qnorm(c(0.25,0.5,0.75)), color="blue")} + 
            geom_boxplot(na.rm=TRUE) + 
            labs(y="sqrt(weight) * residual", x="") +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    }
}





