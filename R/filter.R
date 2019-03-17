
#' Filter a weitrix to ensure each row and column contains sufficient observations
#'
#' When finding components with \code{weitrix_components( )} we will typically need
#'     sufficient observations in each row and column.
#' When finding p components, 
#'     requiring n=p+1 is an absolute minimum, 
#'     and a larger n may reduce the number of random outliers.
#'
#' @param weitrix A weitrix, or something that can be cast to a weitrix with \code{as_weitrix( )}.
#'
#' @param n Number of observations required in each row and in each column.
#'
#' @param weigh_thresh Threshold weight. Can be changed from zero if weights smaller than a certain amount should be disregarded as observations.
#'
#' @export
weitrix_filter <- function(weitrix, n, weight_thresh=0) {
    weitrix <- as_weitrix(weitrix)
    present <- weitrix_weights(weitrix) > weight_thresh

    keep_row <- rep(T, nrow(weitrix))
    keep_col <- rep(T, ncol(weitrix))
    total <- sum(keep_row) + sum(keep_col)
    repeat {
        keep_row[keep_row] <- rowSums(present[keep_row, keep_col, drop=F]) >= n
        keep_col[keep_col] <- colSums(present[keep_row, keep_col, drop=F]) >= n
        new_total <- sum(keep_row) + sum(keep_col)
        if (new_total == total) break
        total <- new_total
    }

    weitrix[keep_row, keep_col]
}

