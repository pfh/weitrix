
#' Weitrix of proportions of totals
#'
#' Construct a weitrix of proportions by dividing a matrix of counts by a matrix of totals. Each proportion is given the total as its weight. The reasoning is that a proportion is an average of "count" 1s and "total-count" 0s. When values are averages, the variance is inversely proportional to the number of individual items being averaged.
#'
#' A reasonable default formula is set in the metadata for use with \code{weitrix_calibrate_trend}. This uses two columns included in the row data. These are the total of totals ("total_reads"), and the expected variance if proportions are uniform ("per_read_var"). 
#'
#' Where the total is zero, the weight will be zero and the proportion will be NA.
#'
#' Counts and totals do not need to be integers, as long as the basic idea of variance being inversely proportional to the "total" holds.
#'
#' @param counts A matrix (or DelayedArray) of counts.
#'
#' @param totals A matrix (or DelayedArray) of totals, of the same dimensions as the counts matrix.
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#'
#' @export
counts_in_totals_proportions <- function(counts, totals) {
    assert_that( identical(dim(counts), dim(totals)) )

    proportions <- counts / totals
    proportions[totals == 0] <- NA
    proportions <- realize_if_delayed(proportions)

    total_counts <- rowSums(counts)
    total_totals <- rowSums(totals)
    mean_props <- total_counts/total_totals
    per_read_var <- mean_props*(1-mean_props)

    result <- as_weitrix(proportions, weights=totals)
    rowData(result)$total_reads <- total_totals
    rowData(result)$per_read_var <- per_read_var

    if (min(per_read_var)/max(per_read_var) < 0.99999999)
        metadata(result)$weitrix$trend_formula <- 
            "~log(per_read_var)+splines::ns(log(total_reads),3)"
    else
        metadata(result)$weitrix$trend_formula <- 
            "~splines::ns(log(total_reads),3)"

    result 
}




