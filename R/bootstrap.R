
# For procrustes rotation, see
# https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem

# TODO: For designed experiments (as opposed to samples from a population), could use:
# - bootstrap within groups, and of blocks
# - parametric bootstrap (need to estimate variance per row)



#' Experimental! Bootstrap resample columns of a weitrix
#' 
#' Bootstrap resample a weitrix. Where a column is sampled multiple times its weights are increased that many times, rather than producing multiple columns.
#'
#' This is experimental, and may be changed or removed in future.
#' 
#' Note that this will not work with methods that count the residual degrees of freedom, such as estimating the residual standard deviation. The main applicaition is with linear model fitting by weighted least squares, and weitrix_components which is based on this. These will work correctly.
#' 
#' @param weitrix A weitrix object, or an object that can be converted to a weitrix with \code{as_weitrix}.
#'
#' @param groups (Optional.) A vector of length ncol(weitrix). Distinct values are treated as groups which are sampled as a single unit. For example the groups could be distinct individuals in a study with "before" and "after" samples.
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


#' Experimental! Bootstrap resample columns of a weitrix then refit components of variation
#' 
#' A bootstrap version of the weitrix is created, and then components refitted by re-optimizing rows, cols, rows, cols, rows, cols, etc (up to \code{max_iter}).
#'
#' This is experimental, and may be changed or removed in future.
#'
#' @param comp A Components object, as produced by weitrix_components.
#' 
#' @param weitrix The weitrix (or an object that can be converted to a weitrix with \code{as_weitrix}) from which comp was produced.
#'
#' @param max_iter Maximum iterations when refitting components.
#'
#' @param verbose Show messages about the progress of the iterations when refitting.
#'
#' @export
components_bootstrap <- function(comp, weitrix, max_iter=5, verbose=TRUE) {
    weitrix <- as_weitrix(weitrix)
    
    if (length(comp$ind_components) == 0) 
        max_iter = 0
    
    booted_weitrix <- weitrix_bootstrap(weitrix)
    x <- weitrix_x(booted_weitrix)
    weights <- weitrix_weights(booted_weitrix)
    
    retained <- metadata(booted_weitrix)$weitrix$retained
    retained_freqs <- metadata(booted_weitrix)$weitrix$retained_freqs
    
    col <- comp$col[retained,,drop=F]
    design <- col[,comp$ind_design,drop=F]
    row <- fit_all_rows(col, x, weights)
    
    for(i in seq_len(max_iter)) {
        if (verbose)
            cat(i,"\n")
        
        col[,comp$ind_components] <- 
            fit_all_cols(row[,comp$ind_components,drop=F], 
                x, weights, row[,comp$ind_design,drop=F], design)

        row <- fit_all_rows(col, x, weights)
    }
    
    rownames(row) <- rownames(x)
    colnames(row) <- colnames(comp$row)
    rownames(col) <- colnames(x)
    colnames(col) <- colnames(comp$col)
    
    comp$col <- col
    comp$row <- row
    comp$R2 <- NULL
    comp$all_R2s <- NULL
    comp
    
    #p <- length(comp$ind_components)
    #design <- comp$col[retained,comp$ind_design,drop=F]
    #initial <- comp$col[retained,comp$ind_components,drop=F]
    #new_comp <- weitrix_components(
    #    booted_weitrix, p=p, design=design, 
    #    n_restarts=1, max_iter=max_iter, use_varimax=FALSE,
    #    initial=initial, verbose=verbose)
    #
    ## Procrustes alignment of col
    #
    #ind <- comp$ind_components
    #M <- crossprod(comp$col[retained,ind,drop=F], new_comp$col[,ind,drop=F] * retained_freqs)
    ## covariance, weighted according to bootstrap
    #decomp <- svd(M)
    #rot <- tcrossprod(decomp$u, decomp$v)
    #
    #new_comp$col[,ind] <- new_comp$col[,ind,drop=F] %*% t(rot)
    #new_comp$row[,ind] <- new_comp$row[,ind,drop=F] %*% t(rot)
    #
    #new_comp
}







