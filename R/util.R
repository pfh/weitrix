
# Utility/helper functions


# Realize a DelayedArray, but don't infect non-delayed arrays
# Checks Delayedness of first argument by default, but can check a second argument instead
realize_if_delayed <- function(x, check=x) {
    if (is(check, "DelayedArray"))
        x <- realize(x)
    x
}



#
# Blocking helper functions
#
# DelayedArray has various utilties for this, but so far I've found them awkward to use.
#
# We will however respect the DelayedArray auto block size
#

# item_size is how many doubles per item
partitions <- function(n, item_size, max_bytes=getAutoBlockSize(), BPPARAM=NULL, cpu_heavy=FALSE) {
    if (n == 0L) return(list())

    step <- max(1L, as.integer(max_bytes/(8*item_size)))
    step <- as.integer(ceiling(n/ceiling(n/step)))
    
    # Overhead of multiprocessing is too high if cluster not already running
    if (cpu_heavy && !is.null(BPPARAM) && bpisup(BPPARAM))
        step <- min(step, as.integer(ceiling(n/bpnworkers(BPPARAM))))
    
    step <- max(1L, step)
    
    starts <- seq(1L, n, by=step)
    ends <- c(starts[-1L]-1L, n)

    map2(starts, ends, seq)
}

# Evaluate an expression with a running worker pool
#
# Otherwise pool will be started for each individual operation!
#
# Only in the default case of a MulticoreParam,
# also sets blas threads in worker pool!
# This may need further tweaking.
#
with_bp_up <- function(expr) {
    BPPARAM <- getAutoBPPARAM()
    if (!bpisup(BPPARAM)) {
        on.exit(bpstop(BPPARAM))
        bpstart(BPPARAM)
        
        if (is(BPPARAM, "MulticoreParam") && bpnworkers(BPPARAM) > 1)
            bplapply(rep(1,bpnworkers(BPPARAM)), blas_set_num_threads, BPPARAM=BPPARAM)
    }

    expr
}


#' Convert a matrix to long form for ggplotting
#'
#' @export
matrix_long <- function(matrix, row_info=NULL, col_info=NULL, varnames=c("name","col")) {
    matrix <- as.matrix(matrix)
    long <- melt(matrix, varnames=varnames)

    if (!is.null(row_info)) {
        name_col <- varnames[1]
        row_info <- as.data.frame(row_info)
        row_info[[name_col]] <- factor(row_info[[name_col]], rownames(matrix))
        long <- left_join(long, row_info, by=name_col)
    }

    if (!is.null(col_info)) {
        name_col <- varnames[2]
        col_info <- as.data.frame(col_info)
        col_info[[name_col]] <- factor(col_info[[name_col]], colnames(matrix))
        long <- left_join(long, col_info, by=name_col)
    }

    long
}











