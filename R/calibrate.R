
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

    if (is(design, "Components"))
        comp <- design
    else
        comp <- weitrix_components(weitrix, design=design, p=0, verbose=FALSE)
    
    assert_that(nrow(comp$row) == nrow(weitrix))
    assert_that(nrow(comp$col) == ncol(weitrix))

    calc_row_dispersion(
        weitrix_x(weitrix),
        weitrix_weights(weitrix),
        comp$row,
        comp$col)        
}


#' Adjust weights based on given row dispersions
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




#' Adjust weights by fitting a trend to estimated dispersions
#'
#' Dispersions are estimated using \code{weitrix_dispersions}. 
#' A trend line is then fitted to log dispersions using a linear model. 
#' Weitrix weights are calibrated based on this trend line. 
#' Any zero or very-near-zero dispersions are ignored when fitting this model.
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
#' If absent, metadata(weitrix)$weitrix$trend_formula is used.
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
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

    if (is.null(trend_formula))
        trend_formula <- metadata(weitrix)$weitrix$trend_formula
    assert_that(!is.null(trend_formula))
    
    trend_formula <- as.formula(trend_formula)

    rowData(weitrix)$dispersion_before <- weitrix_dispersions(weitrix, design)
    
    trend_formula <- update(trend_formula, log(dispersion_before)~.)
    
    data <- rowData(weitrix)
    tiny <- mean(data$dispersion_before,na.rm=TRUE)*1e-9
    good <- !is.na(data$dispersion_before) & data$dispersion_before > tiny
    data <- data[good,]
    
    fit <- eval(substitute(
        lm(trend_formula, data=data),
        list(trend_formula=trend_formula)
    ))
    
    pred <- exp(
        predict(fit, newdata=rowData(weitrix)) +
        sigma(fit)^2/2)
    
    rowData(weitrix)$dispersion_trend <- pred
    rowData(weitrix)$dispersion_after <- 
        rowData(weitrix)$dispersion_before / pred
    
    metadata(weitrix)$weitrix$trend_fit <- fit
    
    weitrix_calibrate(weitrix, pred)
}




