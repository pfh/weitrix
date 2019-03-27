
#
# Blocking helper functions
#
# DelayedArray has various utilties for this, but so far I've found them awkward to use.
#
# We will however respect the DelayedArray auto block size
#

# item_size is how many doubles per item
partitions <- function(n, item_size, max_bytes=getAutoBlockSize()) {
    if (n == 0L) return(list())

    step <- max(1L, as.integer(max_bytes/(8*item_size)))
    step <- max(1L, as.integer(ceiling(n/ceiling(n/step))))
    starts <- seq(1L, n, by=step)
    ends <- c(starts[-1L]-1L, n)

    map2(starts, ends, seq)
}
