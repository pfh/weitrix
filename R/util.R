
# Utility/helper functions


# Realize a DelayedArray, but don't infect non-delayed arrays.
# Checks Delayedness of first argument by default, 
#     but can check a second argument instead.
realize_if_delayed <- function(x, check=x) {
    if (is(check, "DelayedArray"))
        x <- realize(x)
    x
}


#
# Blocking helper functions
#
# DelayedArray has various utilties for this, 
#     but so far I've found them awkward to use.
#
# We will however respect the DelayedArray auto block size
#
# item_size is how many doubles per item
partitions <- function(
        n, item_size, max_bytes=getAutoBlockSize(), 
        BPPARAM=NULL, cpu_heavy=FALSE) {
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


# Start up a worker pool that stays up until the calling function exits.
#
# Otherwise pool will be started for each individual operation!
#
# Sets BLAS threads to 1 before starting worker pool!
# For a fork-based worker pool, 
#     this will ensure each worker only tries to use one BLAS thread.
# BLAS threads are reset to default after the expression is calculated.
# This may need further tweaking.
#
this_function_bp_up <- function() {
    BPPARAM <- bpparam()

    if (is(BPPARAM, "MulticoreParam") && 
            names(dev.cur()) %in% c("quartz","X11cairo","X11")) {
        message(
            "Disabling forking parallel processing as graphics device is open.")
        BPPARAM <- SerialParam()
        BiocParallel::register(BPPARAM)
    }

    # Make sure DelayedArray will use the default bpparam
    # Clean up on function exit
    old_auto_bpparam <- getAutoBPPARAM()
    do.call(
        on.exit,
        list(
            substitute({
                setAutoBPPARAM(old_auto_bpparam)
            }),
            add=TRUE
        ),
        envir = parent.frame()
    )

    setAutoBPPARAM(BPPARAM)


    # Make sure the workers are running
    # Clean up on exit
    if (!bpisup(BPPARAM)) {
        # Deactive X11 graphics before starting forking multiprocessing.
        # May need to add to list of affected devices.

        # Only use 1 thread for linear algebra
        old_threads <- blas_get_num_procs()
        blas_set_num_threads(1)
        
        do.call(
            on.exit, 
            list(
                substitute({
                    bpstop(BPPARAM)
                    blas_set_num_threads(old_threads)
                }), 
                add=TRUE
            ), 
            envir = parent.frame()
        )
        bpstart(BPPARAM)        
    }
}

# Run an expression with a running worker pool.
with_bp_up <- function(expr) {
    this_function_bp_up()
    expr
}


#' Convert a matrix to long form for ggplotting
#'
#' A convenience function which melts the matrix and then 
#'     joins row and column information.
#'
#' @param matrix 
#' A matrix, or object that can be converted to a matrix.
#'
#' @param row_info 
#' Information about rows of the matrix. 
#' A data frame, or object that can be converted to a data frame.
#'
#' @param col_info 
#' Information about columns of the matrix. 
#' A data frame, or object that can be converted to a data frame.
#'
#' @param varnames 
#' Vector of two column names in the output, 
#'     the first for the row and the second for the column.
#'
#' @return
#' A data frame containing the matrix and associated information in long format.
#'
#' @examples
#' matrix_long(weitrix_x(simwei), rowData(simwei), colData(simwei))
#'
#' @export
matrix_long <- function(
        matrix, row_info=NULL, col_info=NULL, varnames=c("name","col")) {
    matrix <- as.matrix(matrix)

    long <- melt(matrix, varnames=varnames, as.is=TRUE)

    if (!is.null(row_info)) {
        name_col <- varnames[1]
        row_info <- as.data.frame(row_info)
        row_info[[name_col]] <- rownames(matrix)
        long <- left_join(long, row_info, by=name_col)
    }

    if (!is.null(col_info)) {
        name_col <- varnames[2]
        col_info <- as.data.frame(col_info)
        col_info[[name_col]] <- colnames(matrix)
        long <- left_join(long, col_info, by=name_col)
    }

    long[[ varnames[1] ]] <- factor(long[[ varnames[1] ]], rownames(matrix))
    long[[ varnames[2] ]] <- factor(long[[ varnames[2] ]], colnames(matrix))

    long
}


