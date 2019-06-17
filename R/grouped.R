
#
# grouped counts converted to proportions
#
# Similar to shifts, but we end up with the same number of rows as we started
#

# This is a conservative weighting scheme
# It would be possible to do better using a prediction other than the row mean
weighted_proportions <- function(mat, name) {
    n <- nrow(mat)
    m <- ncol(mat)
    totals <- colSums(mat)
    total <- sum(totals)
    good <- totals > 0
    
    props <- matrix(0.5, nrow=n,ncol=m)
    weights <- matrix(0, nrow=n,ncol=m)
    
    props[,good] <- t(t(mat[,good]) / totals[good])

    mean_props <- rowSums(mat) / max(1,total)

    # Avoid infinities: as though at least 1 or n-1 reads in each peak
    clip_prop <- min(0.5,1/max(1,total))
    mean_props <- pmax(clip_prop, mean_props)
    mean_props <- pmin(1-clip_prop, mean_props)

    # Bernoulli variance
    per_read_vars <- mean_props*(1-mean_props)
    
    weights[,good] <- outer(1/per_read_vars, totals[good])
    
    rownames(props) <- rownames(mat)
    rownames(weights) <- rownames(mat)
    
    df <- data.frame(
        group=name,
        per_read_weight=1/per_read_vars,
        stringsAsFactors=FALSE)
    rownames(df) <- rownames(mat)

    list(props=props, weights=weights, df=df)
}


counts_proportions_inner <- function(counts, groups) {
    relevant <- unique(unlist(groups))
    counts <- as.matrix(counts[relevant,,drop=F])
    
    results <- map(seq_along(groups), function(i) {
        weighted_proportions(counts[groups[[i]],,drop=FALSE], names(groups)[i])
    })
    rm(counts)
    
    props <- do.call(rbind, map(results, "props"))
    weights <- do.call(rbind, map(results, "weights"))
    df <- do.call(rbind, map(results, "df"))
    rm(results)
    
    SummarizedExperiment(
        assays=list(
            x=realize(props),
            weights=realize(weights)),
        rowData=df)
}

#' Produce a weitrix of proportions within groups
#'
#' Produce a weitrix of proportions between 0 and 1. The input is read counts
#' at a collection of features in a collection of samples. The
#' features need to be grouped, for example by gene. 
#' The proportions will add to 1 within each group.
#' 
#' @param counts A matrix of read counts. Rows are peaks and columns are samples.
#' 
#' @param grouping A data frame defining the grouping of features. Should have a column "group" naming the group and a column "name" naming the feature (corresponding to \code{rownames(counts)}).
#' 
#' @param biovar Produce weights that allow for biological variation. This puts of soft maximum on the effective number of reads, and stops weights from becoming arbitrarily large for large counts.
#' 
#' @param design For biological variation estimation. A design matrix for a linear model, which could account for batches or experimental design. Leaving this with its default will give a conservative choice of weights, and should be safe.
#' 
#' @param p For biological variation estimation. Attempt to find p further components of variation, using \code{weitrix_components}.
#' 
#' @param verbose If TRUE, output some debugging and progress information.
#'
#' @export
counts_proportions <- function(counts, grouping, biovar=TRUE, design=~1, p=0, verbose=TRUE) {
    assert_that(
        is.data.frame(grouping), 
        "group" %in% colnames(grouping), 
        "name" %in% colnames(grouping))
    groups <- split(grouping$name, grouping$group)

    parts <- partitions(length(groups), ncol(counts)/length(groups)*nrow(counts))
    if (verbose)
        message("Calculating proportions in ",length(parts)," blocks")
    
    result <- lapply(parts, function(part) {
        counts_proportions_inner(counts, groups[part])
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
