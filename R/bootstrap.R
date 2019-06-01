
#' Bootstrap resample columns of a weitrix
#' 
#' Bootstrap resample a weitrix. Where a column is sampled multiple times its weights are increased that many times, rather than producing multiple columns.
#' 
#' Note that this will not work with methods that count the residual degrees of freedom, such as estimating the residual standard deviation. The main applicaition is with linear model fitting by weighted least squares, and weitrix_components which is based on this. These will work correctly.
#' 
#' @param weitrix A weitrix object, or an object that can be converted to a weitrix with \code{as_weitrix}.
#'
#' @export
weitrix_bootstrap <- function(weitrix) {
    weitrix <- as_weitrix(weitrix)
    
    n <- ncol(weitrix)
    
    freqs <- tabulate(sample.int(n,replace=TRUE), n)
    
    retain <- freqs != 0
    retain_freqs <- freqs[retain]
    
    weitrix_out <- weitrix[,retain]
    weitrix_weights(weitrix_out) <- sweep(
        weitrix_weights(weirix_out),
        2, retain_freqs, "*")
    weitrix_out
}


#' Bootstrap resample columns of a weitrix then refit components of variation
#' 
#' @param comp A Components object, as produced by weitrix_components.
#' 
#' @param weitrix The weitrix (or an object that can be converted to a weitrix with \code{as_weitrix}) from which comp was produced.
#' 
#components_bootstrap <- function(comp, weitrix) {
#    weitrix <- as_weitrix(weitrix)
#    
#    ...
#    
#}




