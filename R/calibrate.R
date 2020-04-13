
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
    BPPARAM <- getAutoBPPARAM()

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
        design
    else
        weitrix_components(weitrix, design=design, p=0, verbose=FALSE)
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
    
    assert_that(nrow(comp$row) == nrow(weitrix))
    assert_that(nrow(comp$col) == ncol(weitrix))

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
#' If absent, metadata(weitrix)$weitrix$trend_formula is used.
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
        glm(
            trend_formula, 
            data=data,
            weights=deg_free,
            family=quasi(link="log", variance="mu^2"),
            mustart=rep(mustart, n)
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


#' Adjust weights element-wise by fitting a trend to squared residuals
#'
#' This is a very flexible method of calibrating weights.
#' It should be especially useful if your existing weights account for 
#'     technical variation, 
#'     but there is also biological variation.
#' In this case large weights will tend to be overly optimistic, 
#'     and a non-linear transformation of weights is needed.
#' Residuals are found relative to a fitted model.
#' A trend model is then fitted to the squared residuals using a gamma GLM
#'     with log link function.
#' Weitrix weights are set to be the inverse of the fitted trend. 
#'
#' This function is currently not memory efficient,
#'    it should be fine for bulk experiments but not single cell.
#'
#' \code{trend_formula} may reference any row or column variables,
#'     or special factors \code{row} and \code{col},
#'     or \code{weight} for the existing weights.
#' Keep in mind also that a log link function is used.
#' 
#' Unlike in \code{weitrix_calibrate_trend},
#'     existing weights must be explicitly included in the formula
#'     if they are to be retained (see examples).
#' 
#' Example formulae:
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
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#'
#' \code{metadata(weitrix)$weitrix} will contain the fitted trend model,
#'    and if requested the data frame used to fit the model.
#'
#' @examples
#'
#' simcal <- weitrix_calibrate_all(simwei, ~1, ~log(weight))
#'
#' metadata(simcal)$weitrix$all_fit
#'
#' @export
weitrix_calibrate_all <- function(
        weitrix, design=~1, trend_formula=~log(weight), keep_data=FALSE) {

    weitrix <- as_weitrix(weitrix)
    comp <- as_components(design, weitrix)

    # data is very big. Not currently viable for large datasets.
    data <- cbind(
        as.data.frame(rowData(weitrix))[rep(seq_len(nrow(weitrix)), ncol(weitrix)),,drop=FALSE],
        as.data.frame(colData(weitrix))[rep(seq_len(ncol(weitrix)), each=nrow(weitrix)),,drop=FALSE]
    )
    data$row <- 
        rep(factor(rownames(weitrix), rownames(weitrix)), ncol(weitrix))
    data$col <- 
        rep(factor(colnames(weitrix), colnames(weitrix)), each=nrow(weitrix))
    data$weight <- as.vector(weitrix_weights(weitrix))
    data$mu <- as.vector(comp$row %*% t(comp$col))
    data$.y <- (as.vector(weitrix_x(weitrix)) - data$mu)^2
    
    # Want to be able to use poly(log(weight),2)
    # so zero weights need to be dropped.
    good <- data$weight > 0 
    good_data <- data[good,,drop=FALSE]

    mustart <- mean(data$.y, na.rm=TRUE)
    n <- nrow(good_data)

    trend_formula <- update(trend_formula, .y ~ .)
    
    fit <- eval(substitute(
        glm(
            trend_formula,
            data=good_data,
            family=quasi(link="log", variance="mu^2"),
            mustart=rep(mustart, n),
            model=FALSE, x=FALSE, y=FALSE
        ),
        list(trend_formula=trend_formula, mustart=mustart, n=n)
    ))

    pred <- predict(fit, newdata=good_data, type="response")
    pred[is.na(pred)] <- Inf

    weitrix_weights(weitrix)[good] <- 1/pred
    metadata(weitrix)$weitrix$all_fit <- fit

    if (keep_data) {
        data$new_weight <- as.vector(weitrix_weights(weitrix))
        metadata(weitrix)$weitrix$all_data <- data
    }

    weitrix
}





