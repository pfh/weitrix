
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
