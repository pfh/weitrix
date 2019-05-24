
# Shift statistic and associated weights

#
# Converts a matrix of read counts into a vector of shifts relative to the average
# Also provides appropriate weights (1/variance) for technical variation
#
# Samples with less than min_read reads are given shift=0, weight=0
#
weighted_shift <- function(mat, min_reads=1) {
    n <- nrow(mat)
    totals <- colSums(mat)
    good <- totals >= min_reads
    good_mat <- mat[,good,drop=F]
    good_totals <- totals[good]
    
    props <- t(t(good_mat)/good_totals)
    
    mid <- rowMeans(props)
    cummid <- cumsum(mid[-n])
    befores <- c(0,cummid)
    afters <- c(1-cummid, 0)
    pos_score <- befores-afters
    
    # Note: sum(pos_score * mid) == 0, sum(mid) == 1
    per_read_var <- sum(pos_score^2 * mid)

    shifts <- rep(NA_real_, length(good))
    shifts[good] <- colSums(props * pos_score)
    
    weights <- rep(0, length(good))
    weights[good] <- good_totals/per_read_var
    
    list(shifts=shifts, weights=weights, totals=totals, per_read_var=per_read_var)
}


counts_shift_inner <- function(counts, groups, min_reads) {
    relevant <- unique(unlist(groups))
    counts <- as.matrix(counts[relevant,,drop=F])
    
    results <- map(groups, function(members) {
        weighted_shift(counts[members,,drop=FALSE], min_reads=min_reads)
    })
    rm(counts)
    
    shifts <- do.call(rbind, map(results, "shifts"))
    weights <- do.call(rbind, map(results, "weights"))
    #totals <- do.call(rbind, map(results, "totals"))
    per_read_weights <- 1/map_dbl(results, "per_read_var")
    rm(results)

    rownames(shifts) <- names(groups)
    rownames(weights) <- names(groups)
    #rownames(totals) <- names(groups)
    #colnames(shifts) <- colnames(counts)
    #colnames(weights) <- colnames(counts)
    #colnames(totals) <- colnames(counts)
    
    # Drop genes with zero count counts only in a single peak
    good <- is.finite(per_read_weights)
    shifts <- shifts[good,,drop=F]
    weights <- weights[good,,drop=F]
    #totals <- totals[good,,drop=F]
    per_read_weights <- per_read_weights[good]

    SummarizedExperiment(
        assays=list(
            x=realize(shifts),
            weights=realize(weights)),
        rowData=data.frame(
            per_read_weight = per_read_weights))
}


#' Produce a weitrix of shift scores
#'
#' Produce a weitrix of shift scores between -1 and 1. The input is read counts
#' at a collection of peaks (or other features) in a collection of samples. The
#' peaks can be grouped by gene, and are ordered within each gene.
#' 
#' For a particular gene, a shift score measures measures the tendency of reads to be upstrand (negative) or downstrand (positive) of the average over all samples. Shift scores range between -1 and 1. 
#' 
#' @param counts A matrix of read counts. Rows are peaks and columns are samples.
#' 
#' @param grouping A data frame defining the grouping of peaks into genes. Should have a column "group" naming the gene and a column "name" naming the peak (corresponding to \code{rownames(counts)}). Within each group, peak names should be ordered from 5' to 3' position.
#' 
#' @param min_reads Minimum reads to produce a shift score. Where there are fewer than this many reads for a combination of gene and sample, a weight of zero is given.
#' 
#' @param bio_var Produce weights that allow for biological variation.
#' 
#' @param design For biological variation estimation. A design matrix for a linear model, which could account for batches or experimental design. Leaving this with its default will give a conservative choice of weights, and should be safe.
#' 
#' @param p For biological variation estimation. Attempt to find p further components of variation, using \code{weitrix_components}.
#' 
#' @param verbose If TRUE, output some debugging and progress information.
#'
#' @export
counts_shift <- function(counts, grouping, min_reads=1, biovar=TRUE, design=~0, p=0, verbose=TRUE) {
    assert_that(
        is.data.frame(grouping), 
        "group" %in% colnames(grouping), 
        "name" %in% colnames(grouping))
    groups <- split(grouping$name, grouping$group)

    # Only use groups of 2 or more peaks
    good <- map_int(groups, length) >= 2
    groups <- groups[good]

    parts <- partitions(length(groups), ncol(counts)/length(groups)*nrow(counts))
    if (verbose)
        message("Calculating shifts in ",length(parts)," blocks")

    result <- lapply(parts, function(part) {
        counts_shift_inner(counts, groups[part], min_reads)
    })
    result <- do.call(rbind, result)
    colnames(result) <- colnames(counts)
    result <- bless_weitrix(result, "x", "weights")
    
    if (biovar) {
        if (verbose)
            message("Calculating biological noise component")
        result <- weight_balance(result, rowData(result)$per_read_weight, design=design, p=p, verbose=verbose)
        metadata(result)$weitrix$effective_max_reads <- 
            metadata(result)$weitrix$tech_scale / metadata(result)$weitrix$bio_scale
    }
    
    result
}



