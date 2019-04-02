

setClass("Components", representation("list"))

setMethod("show", "Components", function(object) {
    cat("Components are:", paste(colnames(object$row), collapse=", "),"\n")
    cat("$row : ", nrow(object$row),"x",ncol(object$row),"matrix\n")
    cat("$col : ", nrow(object$col),"x",ncol(object$col),"matrix\n")
    cat("$R2  : ", object$R2)
})


# Ensure largest components are first, positive skew on columns of $row
components_order_and_flip <- function(comp) {
    # Ensure largest factor first
    scaling <- sqrt(colMeans(comp$col[,comp$ind_factors,drop=F]^2))
    reordering <- comp$ind_factors[order(scaling, decreasing=TRUE)]
    comp$row[,comp$ind_factors] <- comp$row[,reordering,drop=F]
    comp$col[,comp$ind_factors] <- comp$col[,reordering,drop=F]

    # Ensure positive skew (outliers positive)
    flips <- comp$ind_factors[colSums(comp$row[,comp$ind_factors,drop=F]^3) < 0]
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
    if (length(comp$ind_factors) < 2)
        return(comp)

    rotation <- varimax(comp$row[,comp$ind_factors,drop=FALSE], normalize=FALSE)
    comp$row[,comp$ind_factors] <- 
        comp$row[,comp$ind_factors] %*% rotation$rotmat
    comp$col[,comp$ind_factors] <- 
        comp$col[,comp$ind_factors] %*% rotation$rotmat

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
        present <- wi != 0
        fixed <- as.vector(fixed_row[present,,drop=F] %*% fixed_col[i,])
        least_squares(x[present,,drop=F],wi[present],y[present,i] - fixed)
    })
    
    do.call(rbind, result)
})

fit_all_cols <- function(x,y,w, fixed_row,fixed_col, BPPARAM) {
    if (ncol(x) == 0)
        return(matrix(0, nrow=ncol(y), ncol=0))

    parts <- partitions(ncol(y), nrow(y)*2)
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
        present <- wi != 0
        least_squares(x[present,,drop=F],wi[present],y[i,present])
    })
    
    do.call(rbind, result)
})

fit_all_rows <- function(x,y,w, BPPARAM) {
    if (ncol(x) == 0)
        return(matrix(0, nrow=nrow(y), ncol=0))

    parts <- partitions(nrow(y), ncol(y)*2)
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
        present <- wi != 0
        errors <- as.vector(x[present,i]) - as.vector(row[present,,drop=F] %*% col[i,]) 
        total <- total + sum(errors*errors*wi[present])
    }
    total
})

calc_weighted_ss <- function(x, w, row, col, BPPARAM) {
    parts <- partitions(ncol(x), nrow(x)*2)
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


#' Principal components of a weitrix
#'
#' @export
weitrix_components <- function(
        weitrix, p=2, design=~1, max_iter=100, tol=1e-5, 
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
    ind_factors <- p_design+seq_len(p)

    if (is.null(colnames(design)))
        colnames(design) <- paste0("design", seq_len(p_design))

    col_mat <- design
    if (!is.null(initial))
        col_mat <- cbind(col_mat, initial)
    col_mat <- cbind(col_mat, 
        matrix(rnorm(m*(p_total-ncol(col_mat))), nrow=m))

    # Centering to compare against
    row_mat_center <- fit_all_rows(design, x, weights, BPPARAM=BPPARAM)
    #centered <- x - row_mat_center %*% t(design)
    #ss_total <- sum(centered^2*weights, na.rm=TRUE)
    ss_total <- calc_weighted_ss(x, weights, row_mat_center, design, BPPARAM=BPPARAM)
    R2 <- 0

    if (p == 0) 
        max_iter <- 1

    for(i in seq_len(max_iter)) {
        start <- proc.time()["elapsed"]

        #gc()

        # Esure col_mat[,ind_factors] are orthogonal col_mat[,ind_design]
        col_mat[,ind_factors] <- qr.Q(qr(col_mat))[,ind_factors,drop=F]

        # Update row_mat
        row_mat <- fit_all_rows(col_mat, x, weights, BPPARAM=BPPARAM)

        # Update col_mat
        #centered <- x - row_mat[,seq_len(p_design),drop=F] %*% t(design)
        col_mat <- fit_all_cols(row_mat[,p_design+seq_len(p),drop=F], 
            x, weights, row_mat[,ind_design,drop=F], design, BPPARAM=BPPARAM)
        col_mat <- cbind(design, col_mat)

        # Make decomposition matrices orthogonal
        decomp <- orthogonalize_decomp(row_mat[,ind_factors,drop=F], col_mat[,ind_factors,drop=F])
        row_mat[,ind_factors] <- decomp$rows
        col_mat[,ind_factors] <- decomp$cols

        # Check R^2
        #resid <- x - row_mat %*% t(col_mat)
        #ss_resid <- sum(resid^2*weights, na.rm=TRUE) 
        ss_resid <- calc_weighted_ss(x, weights, row_mat, col_mat, BPPARAM=BPPARAM)
        ratio <- ss_resid/ss_total
        last_R2 <- R2
        R2 <- 1-ratio

        end <- proc.time()["elapsed"]
        if (verbose) {
            message(sprintf("Iter %3d R^2=%7.5f %.1fsec",i,R2,end-start))
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

    result <- new("Components", list(
        row=row_mat, col=col_mat, 
        ind_design=ind_design, 
        ind_factors=ind_factors, 
        R2=R2, 
        df_total=df_total,
        df_model=df_model,
        df_residual=df_residual, 
        iters=i))

    if (use_varimax)
        result <- components_varimax(result)
    else
        result <- components_order_and_flip(result)

    result
}


#' Progressively produce decompositions with varying number of components
#'
#' @export
weitrix_components_seq <- function(
        weitrix, p=10, design=~1, max_iter=100, tol=1e-5, 
        use_varimax=TRUE, verbose=TRUE,
        BPPARAM=getAutoBPPARAM()) {
    weitrix <- as_weitrix(weitrix)
    assert_that(is.number(p))
    assert_that(p >= 1)

    result <- rep(list(NULL), p)

    if (verbose)
        message("Finding ",p," components")

    result[[p]] <- weitrix_components(weitrix, p=p, 
        design=design, max_iter=max_iter, tol=tol, 
        verbose=verbose, BPPARAM=BPPARAM)
    
    for(i in rev(seq_len(p-1))) {
        if (verbose)
            message("Finding ",i," components")

        result[[i]] <- weitrix_components(weitrix, p=i, 
            design=design, max_iter=max_iter, tol=tol, 
            use_varimax=FALSE, verbose=verbose, BPPARAM=BPPARAM,
            initial=result[[i+1]]$col[, result[[i+1]]$ind_factors[seq_len(i)]])
    }

    if (use_varimax)
        result <- map(result, components_varimax)

    result
}


#' Proportion more variance explained by adding components one at a time
#'
#' @export
components_seq_scree <- function(comp_seq) {
    R2 <- map_dbl(comp_seq, "R2")
    R2 - c(0, R2[-length(R2)])
}






