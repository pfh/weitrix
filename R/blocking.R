
#
# Blocking helper functions
#
# DelayedArray has various utilties for this, but so far I've found them awkward to use.
#
# We will however respect the DelayedArray auto block size
#

# item_size is how many doubles per item
partitions <- function(n, item_size, max_bytes=getAutoBlockSize(), BPPARAM=NULL) {
    if (n == 0L) return(list())

    step <- max(1L, as.integer(max_bytes/(8*item_size)))
    step <- as.integer(ceiling(n/ceiling(n/step)))
    
    # Overhead of multiprocessing is high
    #if (!is.null(BPPARAM))
    #    step <- min(step, as.integer(ceiling(n/bpnworkers(BPPARAM))))
    step <- max(1L, step)
    
    starts <- seq(1L, n, by=step)
    ends <- c(starts[-1L]-1L, n)

    map2(starts, ends, seq)
}
