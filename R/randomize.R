

weitrix_randomize_inner <- function(w) {
    result <- matrix(0, nrow=nrow(w),ncol=ncol(w))
    present <- which(as.vector(w > 0))
    result[present] <- rnorm(length(present)) / sqrt(w[present])
    realize(result)
}

#' Generate a random normally distributed version of a weitrix
#'
#' Weights are used to choose the standard deviation of values generated.
#'
#' @export
weitrix_randomize <- function(weitrix, BPPARAM=getAutoBPPARAM()) {
    weitrix <- as_weitrix(weitrix)

    #weitrix_x(weitrix) <- 
    #    matrix(rnorm(nrow(weitrix)*ncol(weitrix)), nrow=nrow(weitrix)) /
    #    sqrt(weitrix_weights(weitrix))
    
    w <- weitrix_weights(weitrix)
    parts <- partitions(ncol(w), nrow(w))
    feed <- map(parts, function(part) {
        w[,part,drop=F]
    })
    result <- bplapply(feed, weitrix_randomize_inner, BPPARAM=BPPARAM)
    weitrix_x(weitrix) <- do.call(rbind, result)

    weitrix
}

