
hill_1 <- function(x) {
    x <- x[x>0]
    if (length(x) == 0) return(0)
    p <- x/sum(x)
    exp(-sum(p*log(p)))
}


#' Calculate Hill numbers (effective number of observations) for rows or columns
#'
#' Effective numbers of observations. 
#' \code{order=0} produces count of non-zero weights. 
#' \code{order=1} produces exp(entropy). 
#' \code{order=2} produces the inverse Simpson index.
#'
#' @param weitrix A weitrix object.
#'
#' @param what Calculate for rows ("row") (default) or columns ("col")?
#'
#' @param order Order of the Hill numbers.
#'
#' @return
#' A numeric vector of effective numbers of observations.
#'
#' @examples
#' weitrix_weights(simwei)
#'
#' weitrix_hill(simwei, what="row", order=0)
#' weitrix_hill(simwei, what="row", order=1)
#' weitrix_hill(simwei, what="row", order=2)
#'
#' weitrix_hill(simwei, what="col", order=0)
#' weitrix_hill(simwei, what="col", order=1)
#' weitrix_hill(simwei, what="col", order=2)
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

