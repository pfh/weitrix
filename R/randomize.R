
#' Generate a random normally distributed version of a weitrix.
#'
#' Weights are used to choose the standard deviation of values generated.
#'
#' @export
weitrix_randomize <- function(weitrix) {
    weitrix <- as_weitrix(weitrix)

    weitrix_x(weitrix) <- 
        matrix(rnorm(nrow(weitrix)*ncol(weitrix)), nrow=nrow(weitrix)) /
        sqrt(weitrix_weights(weitrix))
    weitrix
}

