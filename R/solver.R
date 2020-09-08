
# Functions to do with (contrained) least squares



# Least squares, always producing an answer
# Returns a function, which will cache last weighting used for future calls
# Where multiple solutions exist, one should be chosen arbitrarily
# NA is treated as 0 (should be given weight 0 or not present in any case!)
least_squares_func <- function(A, lower, upper) {
    # Evaluate arguments
    A
    lower
    upper

    # ==== Constrained ====
    if (lower > -Inf || upper < Inf) {
        if (!require(osqp, quietly=TRUE))
            stop("Package osqp is required for constrained least squares.")
        
        l <- rep(lower, nrow(A))
        u <- rep(upper, nrow(A))
        pars <- osqp::osqpSettings(
            eps_abs=1e-9, eps_rel=1e-9, 
            verbose=FALSE)

        func <- function(present,w,b, offset) {
            # Scaling of weights doesn't matter, normalize out
            maxw <- max(w)
            if (maxw > 0)
                w <- w/maxw

            b[is.na(b)] <- 0.0
            sw <- sqrt(w)
            swA <- A[present,,drop=FALSE] * sw
            swb <- sw * (b - offset[present]) 
            
            P <- crossprod(swA,swA)
            q <- -crossprod(swA,swb)
            result <- osqp::solve_osqp(P,q,A,l-offset,u-offset, pars=pars)$x
            
            pred <- A %*% result + offset
            if (any(pred < 0.001/2)) {
                print(pred)
                browser()
            }

            result
        }

        return(func)
    }

    # ==== Unconstrained ====
    state <- new.env()
    state$last_present <- NULL
    state$last_w <- NULL
    state$solver <- NULL

    function(present,w,b, offset) {
        # No data or no non-zero weights?
        if (max(0,w) == 0) 
            return(rep(0, ncol(A)))

        # Scaling of weights doesn't matter, normalize out
        w <- w/max(w)

        # Treat as equal if within 1e-9
        if (!identical(state$last_present, present) ||
            max(abs(state$last_w-w)) >= 1e-9) {
            state$last_present <- present
            state$last_w <- w

            # Soldier on:
            # - Missing values become NA
            # - (near)zero singular values nuked            
            sw <- sqrt(w)
            decomp <- svd(A[present,,drop=FALSE]*sw)
            good <- decomp$d > 1e-9*max(decomp$d)
            state$solver <- 
                decomp$v[,good,drop=FALSE] %*% 
                (t(decomp$u[,good,drop=FALSE]*sw)/decomp$d[good])
        }

        b <- b - offset[present]
        b[is.na(b)] <- 0.0
        as.vector(state$solver %*% b)
    }
}


fit_all_cols_inner <- function(args) {
    y <- as.matrix(args$y)
    w <- as.matrix(args$w)
    solver <- args$least_squares_func(args$x, args$lower, args$upper)

    result <- lapply(seq_len(ncol(y)), function(i) {
        wi <- w[,i]
        present <- wi > 0
        fixed <- as.vector(args$fixed_row %*% args$fixed_col[i,])
        solver(present,wi[present],y[present,i], offset=fixed)
    })
    
    do.call(rbind, result)
}

fit_all_cols <- function(x,y,w, fixed_row,fixed_col, lower, upper) {
    if (ncol(x) == 0)
        return(matrix(0, nrow=ncol(y), ncol=0))
        
    BPPARAM <- bpparam()

    parts <- partitions(ncol(y), nrow(y)*2, BPPARAM=BPPARAM, cpu_heavy=TRUE)
    #cat("cols",length(parts),"\n")
    feed <- map(parts, function(part) {
        list(
            x=x,
            y=y[,part,drop=FALSE],
            w=w[,part,drop=FALSE],
            fixed_row=fixed_row,
            fixed_col=fixed_col[part,,drop=FALSE],
            least_squares_func=least_squares_func,
            lower=lower,
            upper=upper)
    })

    result <- bplapply(feed, fit_all_cols_inner, BPPARAM=BPPARAM)

    do.call(rbind, result)
}


fit_all_rows_inner <- function(args) {
    y <- as.matrix(args$y)
    w <- as.matrix(args$w)
    solver <- args$least_squares_func(args$x, args$lower, args$upper)
    offset <- rep(0, nrow(args$x))

    result <- lapply(seq_len(nrow(y)), function(i) {
        wi <- w[i,]
        present <- wi > 0
        solver(present,wi[present],y[i,present], offset=offset)
    })
    
    do.call(rbind, result)
}

fit_all_rows <- function(x,y,w, lower,upper) {
    if (ncol(x) == 0)
        return(matrix(0, nrow=nrow(y), ncol=0))

    BPPARAM <- bpparam()

    parts <- partitions(nrow(y), ncol(y)*2, BPPARAM=BPPARAM, cpu_heavy=TRUE)
    feed <- map(parts, function(part) {
        list(
            x=x,
            y=y[part,,drop=FALSE],
            w=w[part,,drop=FALSE],
            least_squares_func=least_squares_func,
            lower=lower,
            upper=upper)
    })
    #cat("rows",length(parts),"\n")
    result <- bplapply(feed, fit_all_rows_inner, BPPARAM=BPPARAM)

    do.call(rbind, result)
}
