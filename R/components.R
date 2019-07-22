

setClass("Components", representation("list"))

setMethod("show", "Components", function(object) {
    cat("Components are:", paste(colnames(object$row), collapse=", "),"\n")
    cat("$row : ", nrow(object$row),"x",ncol(object$row),"matrix\n")
    cat("$col : ", nrow(object$col),"x",ncol(object$col),"matrix\n")
    cat("$R2  : ", object$R2,"\n")
})


# Ensure largest components are first, positive skew on columns of $row
components_order_and_flip <- function(comp) {
    # Ensure largest factor first
    scaling <- sqrt(colMeans(comp$col[,comp$ind_components,drop=F]^2))
    reordering <- comp$ind_components[order(scaling, decreasing=TRUE)]
    comp$row[,comp$ind_components] <- comp$row[,reordering,drop=F]
    comp$col[,comp$ind_components] <- comp$col[,reordering,drop=F]

    # Ensure positive skew (outliers positive)
    flips <- comp$ind_components[colSums(comp$row[,comp$ind_components,drop=F]^3) < 0]
    comp$row[,flips] <- -comp$row[,flips,drop=F]
    comp$col[,flips] <- -comp$col[,flips,drop=F]

    comp
}

#' Seek meaningful components by varimax rotation
#'
#' Varimax rotation rotates components in a decomposition so that
#'     each component has a few large loadings and the rest small.
#' By the grace of sparsity we hope that these components 
#'     are then individually interpretable.
#'
#' @param comp A Components object.
#'
#' @export
components_varimax <- function(comp) {
    if (length(comp$ind_components) < 2)
        return(comp)

    rotation <- varimax(comp$row[,comp$ind_components,drop=FALSE], normalize=FALSE)
    comp$row[,comp$ind_components] <- 
        comp$row[,comp$ind_components] %*% rotation$rotmat
    comp$col[,comp$ind_components] <- 
        comp$col[,comp$ind_components] %*% rotation$rotmat

    components_order_and_flip(comp)
}


scale_cols <- function(A,s) {
    assert_that(ncol(A) == length(s))
    t(t(A)*s)
}

# Least squares
# Where multiple solutions exist, one should be chosen arbitrarily
least_squares <- function(A,w,b) {
    if (length(b) == 0) 
        return(rep(0, ncol(A)))

    sw <- sqrt(w)
    # qr.solve(A*sw,b*sw)

    decomp <- svd(A*sw)
    as.vector( decomp$v %*% ((t(decomp$u) %*% (b*sw))/decomp$d) )

    ## May be slightly faster, 
    ## but fails rather than choosing an arbitrary minimum
    # WA <- A*w
    # good <- apply(WA != 0, 2, any)
    # WA <- WA[,good,drop=F]
    # A <- A[,good,drop=F]
    # tAWA <- crossprod(WA, A)
    # tAWb <- as.vector(crossprod(WA,b))
    # L <- chol(tAWA)

    # result <- rep(0, length(good))
    # result[good] <- backsolve(L, forwardsolve(L, tAWb, upper.tri=TRUE, transpose=TRUE))
    # result
}

fit_all_cols_inner <- function(args) with(args, {
    y <- as.matrix(y)
    w <- as.matrix(w)

    result <- lapply(seq_len(ncol(y)), function(i) {
        wi <- w[,i]
        present <- wi > 0
        fixed <- as.vector(fixed_row[present,,drop=F] %*% fixed_col[i,])
        least_squares(x[present,,drop=F],wi[present],y[present,i] - fixed)
    })
    
    do.call(rbind, result)
})

fit_all_cols <- function(x,y,w, fixed_row,fixed_col, BPPARAM) {
    if (ncol(x) == 0)
        return(matrix(0, nrow=ncol(y), ncol=0))

    parts <- partitions(ncol(y), nrow(y)*2, BPPARAM=BPPARAM, cpu_heavy=TRUE)
    #cat("cols",length(parts),"\n")
    feed <- map(parts, function(part) {
        list(
            x=x,
            y=y[,part,drop=F],
            w=w[,part,drop=F],
            fixed_row=fixed_row,
            fixed_col=fixed_col[part,,drop=F],
            least_squares=least_squares)
    })

    result <- bplapply(feed, fit_all_cols_inner, BPPARAM=BPPARAM)

    do.call(rbind, result)
}

fit_all_rows_inner <- function(args) with(args, {
    y <- as.matrix(y)
    w <- as.matrix(w)

    result <- lapply(seq_len(nrow(y)), function(i) {
        wi <- w[i,]
        present <- wi > 0
        least_squares(x[present,,drop=F],wi[present],y[i,present])
    })
    
    do.call(rbind, result)
})

fit_all_rows <- function(x,y,w, BPPARAM) {
    if (ncol(x) == 0)
        return(matrix(0, nrow=nrow(y), ncol=0))

    parts <- partitions(nrow(y), ncol(y)*2, BPPARAM=BPPARAM, cpu_heavy=TRUE)
    feed <- map(parts, function(part) {
        list(
            x=x,
            y=y[part,,drop=F],
            w=w[part,,drop=F],
            least_squares=least_squares)
    })
    #cat("rows",length(parts),"\n")
    result <- bplapply(feed, fit_all_rows_inner, BPPARAM=BPPARAM)

    do.call(rbind, result)
}

# Turn decomposition rows %*% t(cols) into an orthogonal version
# Furthermore, the columns of rows will have unit variance
orthogonalize_decomp <- function(rows,cols) {
    if (ncol(rows) == 0)
        return(list(rows=rows,cols=cols))

    n <- nrow(rows)
    decomp <- svd(rows)
    rows <- decomp$u * sqrt(n)
    cols <- scale_cols(cols %*% decomp$v, decomp$d / sqrt(n))
    decomp <- svd(cols)
    cols <- scale_cols(decomp$u, decomp$d)
    rows <- rows %*% decomp$v
    list(rows=rows,cols=cols)
}



# Weighted sum of squares error of a matrix decomposition
calc_weighted_ss_inner <- function(args) with(args, {
    x <- as.matrix(x)
    w <- as.matrix(w)

    total <- 0
    for(i in seq_len(ncol(x))) {
        wi <- w[,i]
        present <- wi > 0
        errors <- as.vector(x[present,i]) - as.vector(row[present,,drop=F] %*% col[i,]) 
        total <- total + sum(errors*errors*wi[present])
    }
    total
})

calc_weighted_ss <- function(x, w, row, col, BPPARAM) {
    parts <- partitions(ncol(x), nrow(x)*2, BPPARAM=BPPARAM)
    feed <- map(parts, function(part) {
        list(
            x=x[,part,drop=F],
            w=w[,part,drop=F],
            row=row,
            col=col[part,,drop=F])
    })
    #cat("ss",length(parts),"\n")
    result <- bplapply(feed, calc_weighted_ss_inner, BPPARAM=BPPARAM)
    sum(unlist(result))
}

## Test:
# n <- 10
# m <- 7
# p <- 3
# A <- matrix(rnorm(n*p),nrow=n)
# B <- matrix(rnorm(m*p),nrow=m)
# decomp <- orthogonalize_decomp(A,B)
# (A %*% t(B)) - (decomp$rows %*% t(decomp$cols))



weitrix_components_inner <- function(
        weitrix, x, weights, p, p_design, design, max_iter, col_mat, 
        ind_components, ind_design, ss_total, tol, verbose, BPPARAM, outer_iter) {
    R2 <- -Inf

    for(i in seq_len(max_iter)) {
        start <- proc.time()["elapsed"]

        #gc()

        # Esure col_mat[,ind_components] are orthogonal col_mat[,ind_design]
        col_mat[,ind_components] <- qr.Q(qr(col_mat))[,ind_components,drop=F]

        # Update row_mat
        row_mat <- fit_all_rows(col_mat, x, weights, BPPARAM=BPPARAM)

        # Update col_mat
        #centered <- x - row_mat[,seq_len(p_design),drop=F] %*% t(design)
        col_mat <- fit_all_cols(row_mat[,p_design+seq_len(p),drop=F], 
            x, weights, row_mat[,ind_design,drop=F], design, BPPARAM=BPPARAM)
        col_mat <- cbind(design, col_mat)

        # Make decomposition matrices orthogonal
        decomp <- orthogonalize_decomp(row_mat[,ind_components,drop=F], col_mat[,ind_components,drop=F])
        row_mat[,ind_components] <- decomp$rows
        col_mat[,ind_components] <- decomp$cols

        # Check R^2
        #resid <- x - row_mat %*% t(col_mat)
        #ss_resid <- sum(resid^2*weights, na.rm=TRUE) 
        ss_resid <- calc_weighted_ss(x, weights, row_mat, col_mat, BPPARAM=BPPARAM)
        ratio <- ss_resid/ss_total
        last_R2 <- R2
        R2 <- 1-ratio

        end <- proc.time()["elapsed"]
        if (verbose && (i <= 5 || i %% 10 == 0)) {
            message(sprintf("Iter %2d/%3d R^2=%7.5f %.1fsec",outer_iter, i,R2,end-start))
        }

        if (R2 - last_R2 <= tol) 
            break
    }


    # Use original row and column names
    # Give factors in the decomposition meaningful names
    rownames(row_mat) <- rownames(weitrix)
    rownames(col_mat) <- colnames(weitrix)
    colnames(row_mat) <- c(colnames(design), 
        map_chr(seq_len(p), ~paste0("C",.)))
    colnames(col_mat) <- colnames(row_mat)

    new("Components", list(
        row=row_mat, col=col_mat, 
        ind_design=ind_design, 
        ind_components=ind_components, 
        R2=R2, 
        iters=i))
}


#' Principal components of a weitrix
#'
#' Finds principal components of a weitrix.
#' If varimax rotation is enabled, these are then rotated to enhance interpretability.
#'
#' An iterative method is used, starting from a random initial set of components.
#' It is possible for this to get stuck at a local minimum.
#' To ameliorate this, the iteration is run \code{n_restarts} times and the best result used.
#' Examine \code{all_R2s} in the output to see if this is happening -- 
#' if the values are not all nearly identital, the iteration is sometimes getting stuck at local minima..
#' Increase \code{n_restarts} to increase the odds of finding the global minimum.
#'
#' No automatic scaling of the weitrix is performed --
#' all the measurements in the weitrix are assumed to be on a comparable scale.
#'
#' @param weitrix A weitrix object, or an object that can be converted to a weitrix with \code{as_weitrix}.
#' @param p Number of components to find.
#' @param design A formula referring to \code{colData(weitrix)} or a matrix, giving predictors of a linear model for the experimental design. By default only an intercept term is used, ie rows are centered before finding components. A more complex formula might be used to account for batch effects. \code{~0} can be used if rows are already centered.
#' @param n_restarts Number of restarts of the iteration to use.
#' @param max_iter Maximum iterations.
#' @param tol Stop iterating if R-squared increased by less than this amount in an iteration.
#' @param use_varimax Use varimax rotation to enhance interpretability of components.
#' @param initial Optional, an initial guess for column components (scores). Can have fewer columns than \code{p}, in which remaining components are initialized randomly. Only used in the first restart.
#' @param verbose Show messages about the progress of the iterations.
#' @param BPPARAM BiocParallel parameters to use. Defaults to the DelayedArray package's automatic setting.
#' @return
#' A "Components" object with the following elements accessible using \code{$}.
#' \describe{
#'     \item{row}{Row matrix, aka loadings. Rows are rows in the weitrix, and columns contain the experimental design (usually just an intercept term), and components.}
#'     \item{col}{Column matrix, aka scores. Rows are columns in the weitrix, and columns contain fitted coefficients for the experimental design, and components.}
#'     \item{R2}{Weighted R squared statistic. The proportion of total variance explained by the components.}
#'     \item{all_R2s}{R2 statistics from all restarts. This can be used to check how consistently the iteration finds optimal components.}
#'     \item{ind_design}{Column indices associated with experimental design.}
#'     \item{ind_components}{Column indices associated with components.}
#' }
#'
#' For a result \code{comp}, the original measurements are approximated by \code{comp$row %*% t(comp$col)}.
#'
#' \code{weitrix_components_seq} returns a list of Components objects, with increasing numbers of components from 1 up to p.
#'
#' @examples
#' # Variables in rows, observations in columns, as is the Bioconductor convention
#' dat <- t(iris[,1:4])
#'
#' # Find two components
#' comp <- weitrix_components(dat, p=2, max_iter=5, n_restart=1)
#'
#' # Examine row and col matrices
#' pairs(comp$row, panel=function(x,y) text(x,y,rownames(comp$row)))
#' pairs(comp$col)
#'
#' @describeIn weitrix_components
#' Find a matrix decomposition with the specified number of components.
#' @export
weitrix_components <- function(
        weitrix, p, design=~1, n_restarts=1, max_iter=100, tol=1e-5, 
        use_varimax=TRUE, initial=NULL, verbose=TRUE,
        BPPARAM=getAutoBPPARAM()) {
    weitrix <- as_weitrix(weitrix)
    assert_that(is.number(p), p >= 0)

    x <- weitrix_x(weitrix)
    weights <- weitrix_weights(weitrix)

    # TODO: chunking
    #x <- as.matrix(x)
    #weights <- as.matrix(weights)

    if (is(design, "formula"))
        design <- model.matrix(design, data=colData(weitrix))

    rownames(x) <- NULL
    colnames(x) <- NULL
    rownames(weights) <- NULL
    colnames(weights) <- NULL

    n <- nrow(x)
    m <- ncol(x)
    p_design <- ncol(design)
    p_total <- p_design+p

    assert_that(nrow(design) == m, n >= p_total, m >= p)

    df_total <- sum(weights > 0)
    df_model <- n*p_total-m*p
    df_residual <- df_total-df_model

    ind_design <- seq_len(p_design)
    ind_components <- p_design+seq_len(p)

    if (is.null(colnames(design)))
        colnames(design) <- map_chr(seq_len(p_design), ~paste0("design",.))

    # Centering to compare against
    row_mat_center <- fit_all_rows(design, x, weights, BPPARAM=BPPARAM)
    #centered <- x - row_mat_center %*% t(design)
    #ss_total <- sum(centered^2*weights, na.rm=TRUE)
    ss_total <- calc_weighted_ss(x, weights, row_mat_center, design, BPPARAM=BPPARAM)

    if (p == 0) {
        max_iter <- 1
        n_restarts <- 1
    }

    result <- NULL
    best_R2 <- -Inf
    R2s <- numeric(n_restarts)
    for(i in seq_len(n_restarts)) {
        col_mat <- design
        if (!is.null(initial)) {
            # Use first columns of initial on initial run
            # In restarts, use a random projection of all columns from initial
            proj_out <- min(ncol(initial), p)
            if (i == 1)
                proj <- diag(ncol(initial))[,seq_len(proj_out),drop=F]
            else
                proj <- matrix(rnorm(proj_out*ncol(initial)), ncol=proj_out)
            col_mat <- cbind(col_mat, initial %*% proj)
        }
        col_mat <- cbind(col_mat, 
            matrix(rnorm(m*(p_total-ncol(col_mat))), nrow=m))

        this_result <- weitrix_components_inner(
            weitrix=weitrix, x=x, weights=weights, p=p, p_design=p_design, design=design, max_iter=max_iter, col_mat=col_mat, 
            ind_components=ind_components, ind_design=ind_design, ss_total=ss_total, tol=tol, 
            verbose=verbose, BPPARAM=BPPARAM, outer_iter=i)
        R2s[i] <- this_result$R2
        if (this_result$R2 > best_R2) {
            result <- this_result
            best_R2 <- this_result$R2
        }
    }

    result$df_total <- df_total
    result$df_model <- df_model
    result$df_residual <- df_residual
    result$all_R2s <- R2s

    if (use_varimax)
        result <- components_varimax(result)
    else
        result <- components_order_and_flip(result)

    result
}


#' @describeIn weitrix_components
#' Produce a sequence of weitrix decompositions with 1 to p components.
#' @export
weitrix_components_seq <- function(
        weitrix, p, design=~1, n_restarts=1, max_iter=100, tol=1e-5, 
        use_varimax=TRUE, verbose=TRUE,
        BPPARAM=getAutoBPPARAM()) {
    weitrix <- as_weitrix(weitrix)
    assert_that(is.number(p))
    assert_that(p >= 1)

    result <- rep(list(NULL), p)

    if (verbose)
        message("Finding ",p," components")

    result[[p]] <- weitrix_components(weitrix, p=p, 
        design=design, n_restarts=n_restarts, max_iter=max_iter, tol=tol, 
        verbose=verbose, BPPARAM=BPPARAM)
    
    for(i in rev(seq_len(p-1))) {
        if (verbose)
            message("Finding ",i," components")

        result[[i]] <- weitrix_components(weitrix, p=i, 
            design=design, n_restarts=n_restarts, max_iter=max_iter, tol=tol, 
            use_varimax=FALSE, verbose=verbose, BPPARAM=BPPARAM,
            initial=result[[i+1]]$col[,result[[i+1]]$ind_components,drop=F])
    }

    if (use_varimax)
        result <- map(result, components_varimax)

    result
}


inner_scree <- function(comp_seq) {
    
}

#' Proportion more variance explained by adding components one at a time
#'
#' Based on the output of \code{components_seq}, work out how much further variance is explained by adding further components.
#'
#' Some possible threshold levels for including further components are also calculated. 
#"
#' The "Kaiser" threshold is chosen based on the assumption that variance will be spread evenly among all possible components. This is slightly optimistic, as some components will explain more variance simply by chance.
#'
#' The "Parallel analysis" threshold is chosen based on varianced explained by a single component in a randomized weitrix. This will generally be somewhat higher than the "Kaiser" threshold.
#'
#' The "Optimistic" thresholds are chosen starting from the "Parallel Analysis" threshold. We view the Parallel Analysis threshold as indicating random variance is split amongst an effective number of samples, which will be somewhat smaller than the real number of samples. As each component is accepted, the pool of remaining variance is reduced by its contribution, and also the number of effective samples is reduced by one. The threshold is then the size of the remaining variance pool divided by the effective remaining number of samples. This is a somewhat ad-hoc method, but may indicate more components are justified than by criteria based on a flat threshold.
#'
#' @param comp_seq A list of Components objects, as produced by \code{components_seq}.
#' @param rand_comp Optional. A Components object with a single component. This should be based on a randomized version of the weitrix, for example as produced by \code{weitrix_components(weitrix_randomize(my_weitrix), p=1)}.
#'
#' @return
#'
#' \code{components_seq_scree} returns a data frame listing the variance explained by each further component.
#'
#' @export
components_seq_scree <- function(comp_seq, rand_comp=NULL) {
    R2 <- map_dbl(comp_seq, "R2")
    n <- length(R2)
    explained <- R2 - c(0, R2[-n])
    
    result <- data.frame(components=seq_len(n), explained=explained)
    result$kaiser <- rep(1/min(nrow(comp_seq[[1]]$row), nrow(comp_seq[[1]]$col)), n)
    if (!is.null(rand_comp)) {
        assert_that(length(rand_comp$ind_components) == 1)
        pa <- rand_comp$R2
        result$pa <- pa

        decay <- rep(0, n)
        pool <- 1.0
        effective_vars <- 1/pa
        for(i in seq_len(n)) {
            decay[i] <- pool/effective_vars
            pool <- pool - explained[i]
            effective_vars <- effective_vars-1
        }
        result$optimistic <- decay
    }
    
    result
}

#' @rdname components_seq_scree
#'
#' @return
#' 
#' \code{components_seq_screeplot} returns a ggplot2 plot object.
#' 
#' @export
components_seq_screeplot <- function(comp_seq, rand_comp=NULL) {
    df <- components_seq_scree(comp_seq, rand_comp)
    n <- nrow(df)
    
    result <- ggplot(df, aes_string(x="components"))
    
    result <- result +
        geom_line(aes_string(y="kaiser", color='"Kaiser Criterion"'))
    
    if (!is.null(df$pa)) {
        result <- result +
            geom_line(aes_string(y="pa", color='"Parallel Analysis"')) +
            geom_line(aes_string(y="optimistic", color='"Optimistic"'))
        
    }
    
    result <- result + 
        geom_point(aes_string(x="components", y="explained")) +
        geom_hline(yintercept=0) +
        labs(x="Number of components", y="Further variance explained",
             color="Threshold guidance")
    result
}



