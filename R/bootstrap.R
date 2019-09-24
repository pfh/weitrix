
# For procrustes rotation, see
# https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem

# TODO: For designed experiments (as opposed to samples from a population), could use:
# - bootstrap within groups, and of blocks
# - parametric bootstrap (need to estimate variance per row)



#' Bootstrap resample columns of a weitrix
#' 
#' Bootstrap resample a weitrix. Where a column is sampled multiple times its weights are increased that many times, rather than producing multiple columns.
#' 
#' Note that this will not work with methods that count the residual degrees of freedom, such as estimating the residual standard deviation. The main applicaition is with linear model fitting by weighted least squares, and weitrix_components which is based on this. These will work correctly.
#' 
#' @param weitrix A weitrix object, or an object that can be converted to a weitrix with \code{as_weitrix}.
#'
#' @export
weitrix_bootstrap <- function(weitrix, groups=NULL) {
    weitrix <- as_weitrix(weitrix)
    
    n <- ncol(weitrix)
    
    if (is.null(groups))
        groups <- seq_len(n)
    assert_that(length(groups) == n)

    groups <- as.integer(factor(groups))
    n_levels <- max(groups)
    
    level_freqs <- tabulate(sample.int(n_levels,replace=TRUE), n_levels)
    freqs <- level_freqs[groups] 
    
    retain <- freqs != 0
    retain_freqs <- freqs[retain]
    
    weitrix_out <- weitrix[,retain]
    metadata(weitrix_out)$weitrix$retained <- retain
    metadata(weitrix_out)$weitrix$retained_freqs <- retain_freqs
    weitrix_weights(weitrix_out) <- sweep(
        weitrix_weights(weitrix_out),
        2, retain_freqs, "*")
    weitrix_out
}


#' Bootstrap resample columns of a weitrix then refit components of variation
#' 
#' @param comp A Components object, as produced by weitrix_components.
#' 
#' @param weitrix The weitrix (or an object that can be converted to a weitrix with \code{as_weitrix}) from which comp was produced.
#' @param max_iter Maximum iterations when refitting components.
#' @param verbose Show messages about the progress of the iterations when refitting.
#'
#' @export
components_bootstrap <- function(comp, weitrix, max_iter=10, verbose=TRUE) {
    weitrix <- as_weitrix(weitrix)
    
    booted_weitrix <- weitrix_bootstrap(weitrix)
    
    retained <- metadata(booted_weitrix)$weitrix$retained
    retained_freqs <- metadata(booted_weitrix)$weitrix$retained_freqs
    
    p <- length(comp$ind_components)
    design <- comp$col[retained,comp$ind_design,drop=F]
    new_comp <- weitrix_components(
        booted_weitrix, p=p, design=design, 
        max_iter=max_iter, use_varimax=FALSE, 
        verbose=verbose)
    
    # Procrustes alignment of col
    
    ind <- comp$ind_components
    M <- crossprod(comp$col[retained,ind,drop=F], new_comp$col[,ind,drop=F] * retained_freqs)
    # covariance, weighted according to bootstrap
    decomp <- svd(M)
    rot <- tcrossprod(decomp$u, decomp$v)
    
    new_comp$col[,ind] <- new_comp$col[,ind,drop=F] %*% t(rot)
    new_comp$row[,ind] <- new_comp$row[,ind,drop=F] %*% t(rot)
    
    new_comp
}







