
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
            name=I(rownames(shifts)),
            per_read_weight = per_read_weights))
}

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
        result <- weight_balance(result, rowData(result)$per_read_weight, design=design, p=p)
        metadata(result)$weitrix$effective_max_reads <- 
            metadata(result)$weitrix$tech_scale / metadata(result)$weitrix$bio_scale
    }
    
    result
}



