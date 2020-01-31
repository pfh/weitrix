
#
# grouped counts converted to proportions
#
# Similar to shifts, but we end up with the same number of rows as we started
#
# Weights are simply total reads
#
weighted_proportions <- function(mat, name) {
    n <- nrow(mat)
    m <- ncol(mat)
    totals <- colSums(mat)
    total <- sum(totals)
    good <- totals > 0
    
    props <- matrix(NA, nrow=n,ncol=m)
    weights <- matrix(0, nrow=n,ncol=m)
    
    props[,good] <- t(t(mat[,good]) / totals[good])

    mean_props <- rowSums(mat) / max(1,total)

    # Bernoulli variance
    per_read_vars <- mean_props*(1-mean_props)
    
    weights[,good] <- outer(rep(1,n), totals[good])
    
    rownames(props) <- rownames(mat)
    rownames(weights) <- rownames(mat)
    
    df <- data.frame(
        group=name,
        total_reads=total,
        per_read_var=per_read_vars,
        stringsAsFactors=FALSE)
    rownames(df) <- rownames(mat)

    list(props=props, weights=weights, df=df)
}


counts_proportions_inner <- function(counts, groups) {
    relevant <- unique(unlist(groups))
    counts <- as.matrix(counts[relevant,,drop=FALSE])
    
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
            x=realize_if_delayed(props),
            weights=realize_if_delayed(weights)),
        rowData=df)
}

#' Produce a weitrix of proportions within groups
#'
#' Produce a weitrix of proportions between 0 and 1. 
#' The input is read counts at a collection of features 
#'     in a collection of samples. 
#' The features need to be grouped, for example by gene. 
#' The proportions will add to 1 within each group.
#' 
#' @param counts 
#' A matrix of read counts. 
#' Rows are peaks and columns are samples.
#' 
#' @param grouping 
#' A data frame defining the grouping of features. 
#' Should have a column "group" naming the group and 
#'     a column "name" naming the feature 
#'     (corresponding to \code{rownames(counts)}).
#' 
#' @param verbose 
#' If TRUE, output some debugging and progress information.
#'
#' @export
counts_proportions <- function(counts, grouping, verbose=TRUE) {
    assert_that(
        is.data.frame(grouping), 
        "group" %in% colnames(grouping), 
        "name" %in% colnames(grouping))
    groups <- split(grouping$name, grouping$group)

    parts <- 
        partitions(length(groups), ncol(counts)/length(groups)*nrow(counts))
    if (verbose)
        message("Calculating proportions in ",length(parts)," blocks")
    
    result <- lapply(parts, function(part) {
        counts_proportions_inner(counts, groups[part])
    })
    result <- do.call(rbind, result)
    colnames(result) <- colnames(counts)
    result <- bless_weitrix(result, "x", "weights")
    metadata(result)$weitrix$trend_formula <- 
        "~log(per_read_var)+splines::ns(log(total_reads),3)"
    result
}
