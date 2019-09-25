
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
            errors <- as.vector(x[i,present]) - as.vector(col[present,,drop=F] %*% row[i,]) 
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
            x=x[part,,drop=F],
            w=w[part,,drop=F],
            row=row[part,,drop=F],
            col=col)
    })
    result <- bplapply(feed, calc_row_dispersion_inner, BPPARAM=BPPARAM)
    unlist(result)
}


#' Calculate row dispersions
#'
#' Calculate the dispersion of each row. For each observation, this value divided by the weight gives the observation's variance.
#'
#' @param weitrix A weitrix object, or an object that can be converted to a weitrix with \code{as_weitrix}.
#' @param comp A Components, an approximate matrix decomposition of the weitrix x values, for example created with \code{weitrix_components}.
#'
#' @export
weitrix_dispersions <- function(weitrix, comp) {
    weitrix <- as_weitrix(weitrix)
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
#' For large numbers of samples this can be based directly on weitrix_dispersions. For small numbers of samples, when using limma, it should be based on a trend-line fitted to known co-variates of the dispersions. This can be done using \code{weitrix_calibrate_trend}.
#'
#' @param weitrix A weitrix object, or an object that can be converted to a weitrix with \code{as_weitrix}.
#' @param dispersions A dispersion for each row.
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
#' Dispersions are estimated using \code{weitrix_dispersions}. A trend line is then fitted to log dispersions using a linear model. Weitrix weights are calibrated based on this trend line. Any zero or very-near-zero dispersions are ignored when fitting this model.
#'
#' @param weitrix A weitrix object, or an object that can be converted to a weitrix with \code{as_weitrix}.
#' @param comp A Components, an approximate matrix decomposition of the weitrix x values, for example created with \code{weitrix_components}.
#' @param formula A formula specification for predicting log dispersion from columns of rowData(weitrix). If absent, metadata(weitrix)$weitrix$trend_formula is used.
#'
#' @export
weitrix_calibrate_trend <- function(weitrix, comp, formula=NULL) {
    weitrix <- as_weitrix(weitrix)
    
    if (is.null(formula))
        formula <- metadata(weitrix)$weitrix$trend_formula
    assert_that(!is.null(formula))
    
    formula <- as.formula(formula)
    
    rowData(weitrix)$dispersion <- weitrix_dispersions(weitrix, comp)
    
    formula <- update(formula, log(dispersion)~.)
    
    data <- rowData(weitrix)
    tiny <- mean(data$dispersion,na.rm=TRUE)*1e-9
    good <- !is.na(data$dispersion) & data$dispersion > tiny
    data <- data[good,]
    
    fit <- eval(substitute(
        lm(formula, data=data),
        list(formula=formula)
    ))
    
    pred <- exp(predict(fit, newdata=rowData(weitrix)))
    
    rowData(weitrix)$dispersion_trend <- pred
    metadata(weitrix)$weitrix$trend_fit <- fit
    
    weitrix_calibrate(weitrix, pred)
}




