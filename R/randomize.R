

weitrix_randomize_inner <- function(w) {
    result <- matrix(0, nrow=nrow(w),ncol=ncol(w))
    present <- which(as.vector(w > 0))
    result[present] <- rnorm(length(present)) / sqrt(w[present])
    
    # Infect with delayedness of w
    realize_if_delayed(result, w)
}

#' Generate a random normally distributed version of a weitrix
#'
#' Values are generated with variance equal to 1/weight.
#' This can be used to see what R-squared would be achieved 
#'     with purely random data,
#' and therefore an appropriate number of components to use.
#' This is known as Parallel Analysis.
#'
#' @param weitrix 
#' A weitrix object, or an object that can be converted to a weitrix 
#'     with \code{as_weitrix}.
#'
#' @return 
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#' 
#' @seealso
#' \code{\link{components_seq_screeplot}}
#'
#' @examples
#' weitrix_randomize(simwei)
#'
#' @export
weitrix_randomize <- function(weitrix) {
    weitrix <- as_weitrix(weitrix)

    #weitrix_x(weitrix) <- 
    #    matrix(rnorm(nrow(weitrix)*ncol(weitrix)), nrow=nrow(weitrix)) /
    #    sqrt(weitrix_weights(weitrix))
    
    w <- weitrix_weights(weitrix)
    parts <- partitions(ncol(w), nrow(w))
    feed <- map(parts, function(part) {
        w[,part,drop=FALSE]
    })
    #Realization currently not working in workers
    #result <- bplapply(feed, weitrix_randomize_inner, BPPARAM=BPPARAM)
    result <- lapply(feed, weitrix_randomize_inner)
    result <- do.call(cbind, result)
    rownames(result) <- rownames(weitrix)
    colnames(result) <- colnames(weitrix)
    weitrix_x(weitrix) <- result

    weitrix
}

