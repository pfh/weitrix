
hill_1 <- function(x) {
    x <- x[x>0]
    if (length(x) == 0) return(0)
    p <- x/sum(x)
    exp(-sum(p*log(p)))
}


#' Calculate Hill numbers (effective number of observations) for rows or columns
#'
#' Effective numbers of observations. order=0 produces count of non-zero weights. order=1 produces exp(entropy). order=2 produces the inverse Simpson index.
#'
#' @export
weitrix_hill <- function(weitrix, what=c("row","col"), order=2) {
    weitrix <- as_weitrix(weitrix)
    what <- match.arg(what)    
    assert_that(order >= 0)
    
    weights <- weitrix_weights(weitrix)
    if (what == "row") 
        weights <- t(weights)
    
    if (order == 0)
        return(colSums(weights>0))
        
    if (order == 1)
        return(apply(weights,2,hill_1))
    
    result <- (colSums(weights)^order/colSums(weights^order))^(1/(order-1))
    result[is.nan(result)] <- 0
    result
}


#' Downweight rows and cols with low effective numbers of observations
#'
#' This is an alternative to filtering out sparse rows and columns, generally to be used before weitrix_components. For each row, weights are scaled by the row's Hill number of order \code{order} as a fraction of the total possible, raised to the power \code{power}. The same is done for each column.
#'
#' @export
weitrix_downweight_sparse <- function(weitrix, power=2, order=2, rows=TRUE, cols=TRUE) {
    weitrix <- as_weitrix(weitrix)
    weitrix_out <- weitrix
    
    if (rows) {
        mult <- (weitrix_hill(weitrix, "row", order)/nrow(weitrix)) ^ power
        weitrix_weights(weitrix_out) <- 
            sweep(weitrix_weights(weitrix_out), 1, mult, "*")
    }

    if (cols) {
        mult <- (weitrix_hill(weitrix, "col", order)/ncol(weitrix)) ^ power
        weitrix_weights(weitrix_out) <- 
            sweep(weitrix_weights(weitrix_out), 2, mult, "*")
    }
    
    weitrix_out
}