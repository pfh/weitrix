
# Shift statistic and associated weights

#
# Converts a matrix of read counts into a vector of shifts 
# relative to the average.
#
# Weights are simply total read count.
#
# Samples with no reads are given shift=NA, weight=0.
#
weighted_shift <- function(mat) {
    n <- nrow(mat)
    totals <- colSums(mat)
    good <- totals > 0
    good_mat <- mat[,good,drop=FALSE]
    good_totals <- totals[good]
    
    props <- t(t(good_mat)/good_totals)

    row_totals <- rowSums(good_mat)
    mid <- row_totals / sum(row_totals)
    # Alternative would be to calculate proportions first
    #mid <- rowMeans(props)    
    
    cummid <- cumsum(mid[-n])
    befores <- c(0,cummid)
    afters <- c(1-cummid, 0)
    pos_score <- befores-afters
    
    # Note: sum(pos_score * mid) == 0, sum(mid) == 1
    per_read_var <- sum(pos_score^2 * mid)

    shifts <- rep(NA_real_, length(good))
    shifts[good] <- colSums(props * pos_score)
    
    weights <- rep(0, length(good))
    weights[good] <- good_totals

    grand_total <- sum(good_totals)
    
    list(
        shifts=shifts, 
        weights=weights, 
        total=grand_total, 
        per_read_var=per_read_var)
}


counts_shift_inner <- function(counts, groups) {
    relevant <- unique(unlist(groups))
    counts <- as.matrix(counts[relevant,,drop=FALSE])
    
    results <- map(groups, function(members) {
        weighted_shift(counts[members,,drop=FALSE])
    })
    rm(counts)
    
    shifts <- do.call(rbind, map(results, "shifts"))
    weights <- do.call(rbind, map(results, "weights"))
    totals <- do.call(rbind, map(results, "total"))
    per_read_vars <- map_dbl(results, "per_read_var")
    rm(results)

    rownames(shifts) <- names(groups)
    rownames(weights) <- names(groups)
    #rownames(totals) <- names(groups)
    #colnames(shifts) <- colnames(counts)
    #colnames(weights) <- colnames(counts)
    #colnames(totals) <- colnames(counts)
    
    # Drop genes with zero count, counts only in a single peak
    good <- is.finite(per_read_vars) & per_read_vars > 0
    shifts <- shifts[good,,drop=FALSE]
    weights <- weights[good,,drop=FALSE]
    totals <- totals[good,,drop=FALSE]
    per_read_vars <- per_read_vars[good]

    SummarizedExperiment(
        assays=list(
            x=realize_if_delayed(shifts),
            weights=realize_if_delayed(weights)),
        rowData=data.frame(
            per_read_var = per_read_vars,
            total_reads = totals))
}


#' Produce a weitrix of shift scores
#'
#' Produce a weitrix of shift scores between -1 and 1. The input is read counts
#' at a collection of peaks (or other features) in a collection of samples. The
#' peaks can be grouped by gene, and are ordered within each gene.
#' 
#' For a particular gene, a shift score measures measures the tendency of 
#' reads to be upstrand (negative) or downstrand (positive) of 
#' the average over all samples. 
#' Shift scores range between -1 and 1. 
#' 
#' @param counts 
#' A matrix of read counts. Rows are peaks and columns are samples.
#' 
#' @param grouping 
#' A data frame defining the grouping of peaks into genes. 
#' Should have a column "group" naming the gene and 
#'     a column "name" naming the peak 
#'     (corresponding to \code{rownames(counts)}). 
#' Within each group, peak names should be ordered from 5' to 3' position.
#' 
#' @param verbose 
#' If TRUE, output some debugging and progress information.
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#' 
#' @examples
#' grouping <- data.frame(
#'     group=c("A","A","A","B","B"),
#'     name=c("p1","p2","p3","p4","p5"))
#'
#' counts <- rbind(
#'     p1=c(1,2,0),
#'     p2=c(0,1,0),
#'     p3=c(1,0,0),
#'     p4=c(0,0,1),
#'     p5=c(0,2,1))
#'
#' wei <- counts_shift(counts, grouping)
#'
#' weitrix_x(wei)
#' weitrix_weights(wei)
#' rowData(wei)
#'
#' @export
counts_shift <- function(counts, grouping, verbose=TRUE) {
    assert_that(
        is.data.frame(grouping), 
        "group" %in% colnames(grouping), 
        "name" %in% colnames(grouping))
    groups <- split(grouping$name, grouping$group)

    # Only use groups of 2 or more peaks
    good <- map_int(groups, length) >= 2
    groups <- groups[good]

    parts <- partitions(
        length(groups), ncol(counts)/length(groups)*nrow(counts))
    if (verbose)
        message("Calculating shifts in ",length(parts)," blocks")

    result <- lapply(parts, function(part) {
        counts_shift_inner(counts, groups[part])
    })
    result <- do.call(rbind, result)
    colnames(result) <- colnames(counts)
    result <- bless_weitrix(result, "x", "weights")
    metadata(result)$weitrix$calibrate_trend_formula <- 
        "~log(per_read_var)+splines::ns(log(total_reads),4)"
    metadata(result)$weitrix$calibrate_all_formula <- 
        "~log(per_read_var)+splines::ns(log(weight),4)"
    
    result
}



