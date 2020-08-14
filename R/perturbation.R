

calc_rowstats_inner <- function(args) with(args, {
    x <- as.matrix(x)
    w <- as.matrix(w)
    n <- nrow(x)
    p <- ncol(row)

    result <- rep(NA_real_,n)

    wss <- rep(NA_real_,n)
    total_weight <- rep(NA_real_,n)
    n_present <- rep(NA_real_,n)
    df <- rep(NA_real_,n)
    
    for(i in seq_len(n)) {
        wi <- w[i,]
        present <- wi > 0
        wi_present <- wi[present]
        n_present[i] <- length(wi_present)
        df[i] <- n_present[i]-p
        total_weight[i] <- sum(wi_present)
        if (df[i] > 0) {
            errors <- as.vector(x[i,present]) - 
                as.vector(col[present,,drop=FALSE] %*% row[i,]) 
            wss[i] <- sum(errors*errors*wi_present)
        }
    }

    data.frame(
        wss=wss,
        total_weight=total_weight,
        df=df,
        n_present=n_present)
})

calc_rowstats <- function(x, w, row, col) {
    BPPARAM <- bpparam()

    parts <- partitions(nrow(x), ncol(x)*2, BPPARAM=BPPARAM)
    feed <- map(parts, function(part) {
        list(
            x=x[part,,drop=FALSE],
            w=w[part,,drop=FALSE],
            row=row[part,,drop=FALSE],
            col=col)
    })
    result <- bplapply(feed, calc_rowstats_inner, BPPARAM=BPPARAM)
    do.call(rbind, result)
}

#' Find rows with confidently excessive variability in a calibrated weitrix
#'
#' Find rows with confident excess standard deviation beyond what is expected based on the weights of a calibrated weitrix. 
#' This may be used, for example, to find potential marker genes.
#'
#' This is a conversion of the "dispersion" statistic for each row into units that are more readily interpretable, accompanied by confidence bounds with a multiple testing correction.
#'
#' We are looking for further perturbation of observed values beyond what is accounted for by a linear model and, further, beyond what is expected based on the observation weights (assumed to be calibrated and so interpreted as 1/variance). We are seeking to estimate the standard deviation of this further perturbation.
#'
#' The weitrix must have been calibrated for results to make sense.
#'
#' Top confident effect sizes are found using the \code{topconfects} method, based on the model that the observed weighted sum of squared residuals being non-central chi-square distributed.
#'
#' Note that all calculations are based on weighted residuals, with a rescaling to place results on the original scale. When a row has highly variable weights, this is an approximation that is only sensible if the weights are unrelated to the values themselves. 
#' 
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#' with \code{as_weitrix}.
#' @param design 
#' A formula in terms of \code{colData(weitrix} or a design matrix, 
#' which will be fitted to the weitrix on each row. 
#' Can also be a pre-existing Components object, 
#' in which case the existing fits (\code{design$row}) are used.
#' @param fdr False Discovery Rate to control for.
#' @param step Granularity of effect sizes to test.
#'
#' @return
#' A topconfects result. The \code{$table} data frame contains columns:
#' \itemize{
#' \item{effect}{ Estimated excess standard deviation, in the same units as the observations themselves. 0 if the dispersion is less than 1. }
#' \item{confect}{ A lower confidence bound on effect. }
#' \item{typical_obs_err}{ Typical accuracy of each observation. }
#' \item{dispersion}{ Dispersion. Weighted sum of squared residuals divided by residual degrees of freedom. }
#' \item{n_present}{ Number of observations with non-zero weight. }
#' \item{df}{ Degrees of freedom. n minus the number of coefficients in the model. }
#' \item{fdr_zero}{ FDR-adjusted p-value for the null hypothesis that effect is zero. }
#' }
#'
#' Note that \code{dispersion = effect^2/typical_obs_err^2 + 1} for non-zero effect values.
#'
#' @examples
#'
#' # weitrix_sd_confects should only be used with a calibrated weitrix
#' calwei <- weitrix_calibrate_all(simwei, ~1, ~1)
#'
#' weitrix_sd_confects(calwei, ~1)
#'
#' @export
weitrix_sd_confects <- function(
        weitrix, design=~1,  
        fdr=0.05, step=0.001) {
    weitrix <- as_weitrix(weitrix)
    comp <- as_components(design, weitrix)
    
    assert_that(nrow(comp$row) == nrow(weitrix))
    assert_that(nrow(comp$col) == ncol(weitrix))

    # Collect some numbers

    df <- calc_rowstats(
        weitrix_x(weitrix),
        weitrix_weights(weitrix),
        comp$row,
        comp$col)
    
    df$dispersion <- df$wss / df$df
    df$mean_weight <- df$total_weight / df$n_present
    df$typical_obs_err <- sqrt(1/df$mean_weight)
    df$effect <- sqrt(pmax(0, (df$dispersion-1)/df$mean_weight ))

    # Do the topconfects

    pfunc <- function(indices, effect_size) {
        1 - pchisq(
            df$wss[indices],
            df$df[indices],
            ncp=effect_size^2 * df$mean_weight[indices] * df$df[indices])
    }

    result <- nest_confects(
        nrow(df), pfunc, fdr=fdr, step=step, full=TRUE)
    
    # Add further useful columns to the results table

    fdr_zero <- result$table$fdr_zero
    result$table$fdr_zero <- NULL
    
    result$table$effect <- df$effect[result$table$index]

    for(name in c("typical_obs_err","dispersion","n_present","df"))
        result$table[[name]] <- df[result$table$index,name]

    result$table$fdr_zero <- fdr_zero

    result$table$name <- rownames(weitrix)[result$table$index]

    for(name in colnames(rowData(weitrix)))
        if (!name %in% colnames(result$table))
            result$table[[name]] <- rowData(weitrix)[result$table$index,name]

    result$effect_desc <- "excess standard deviation"

    result
}
