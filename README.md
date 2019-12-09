# weitrix

R package for working with matrices with precision weights.

Install with:

```
BiocManager::install("pfh/weitrix")
```

* [GitHub page](https://github.com/pfh/weitrix)
* [Documentation site](http://logarithmic.net/weitrix/)

## Weights

A weight determines the importance of an observation. One way of thinking about a weight is that it is as if this one observation is actually an average over some number of real observations. For example if an observation is an average over some number of reads, the number of reads might be used as the weight. 

The choice of weight is somewhat arbitrary. You can use it simply to tell model fitting (such as that in `weitrix_components`) what to pay most attention to. It's better to pay more attention to more accurately measured observations.

A weight of 0 indicates completely missing data.

Weights can be calibrated so they are one over the variance of a measurement. This may be a necessary step, which in this package we will have happen after model fitting but before statistical testing.

## What is a weitrix?

A "weitrix" is a SummarizedExperiment object (or subclass thereof) with two assays, one containing the actual measurements and the other the associated weights. A "weitrix" metadata entry stores the names of these assays. There are several ways to construct a weitrix:

* `as_weitrix(x, weights)` constructs a weitrix, where `x` is a matrix of measurements and `weights` is a corresponding matrix of weights.

* A SummarizedExperiment can be marked as a weitrix using `bless_weitrix`. This requires specifying the names of the two assays to be used.

* Anything the limma package knows how to work with can be converted to a weitrix using `as_weitrix`. (Most functions in the weitrix package will attempt this conversion automatically.)

The usual SummarizedExperiment accessor functions can be used: `assay` `rowData` `colData` `metadata`

Additionally, the blessed assays be accessed using: `weitrix_x` `weitrix_weights`

## Rows and columns

`weitrix` follows the Bioconductor convention of placing features in rows and units of observation (samples, cells) in columns. This is the transpose of the normal R convention!

## Use with limma or topconfects

A weitrix can be converted to an EList object for use with limma: `weitrix_elist`

The `$col` matrix of a `Components` may be used as a design matrix for differential analysis with limma. **Warning:** This may produce liberal results, because the design matrix is itself uncertain and this isn't taken into account. Use this with caution.

## Big datasets

weitrix can use DelayedArray assays. Functions that produce weitrices will used `DelayedArray` output assays if given `DelayedArray` input assays.

weitrix will attempt to perform calculations blockwise in parallel. weitrix tries to use DelayedArray defaults. Adjust with `DelayedArray::setRealizationBackend`, 
`DelayedArray::setAutoBlockSize`, `DelayedArray::setAutoBPPARAM`.

It is always possible to convert an assay back to a normal R matrix with `as.matrix`.

Set the DelayedArray realization backend to `RleMatrix` or `HDF5Array` if weitrices will be too big to fit in memory uncompressed. The `RleMatrix` stores data in memory in a compressed form. The `HDF5Array` backend stores data on disk, in temporary files.

* `DelayedArray::setRealizationBackend("RleArray")` may require copying the complete matrix to all members of the cluster (not sure on this point). HDF5Array seems better, since all the actual data is stored on disk and the disk cache will be shared between processes.

* If using `DelayedArray::setRealizationBackend("HDF5Array")` you may also want to set `HDF5Array::setHDF5DumpDir`.

* A weitrix can be permanently stored to disk using `HDF5Array::saveHDF5SummarizedExperiment`.

Example setup:

```
library(DelayedArray)
library(HDF5Array)

# Store intermediate results in a directory called __dump__
# You may need to clean up this directory manually
setRealizationBackend("HDF5Array")
setHDF5DumpDir("__dump__")
```

## Parallelism fine tuning

### BiocParallel problems

Parallel processing in R and Bioconductor remains finicky but is also necessary for large datasets. weitrix uses DelayedArray's default parallel processing as its own default, see `DelayedArray::getAutoBPPARAM()` and `DelayedArray::setAutoBPPARAM()`.

If weitrix hangs, try running it with serial processing:

```
DelayedArray::setAutoBPPARAM( BiocParallel::SerialParam() )
```

### OpenBLAS

If using parallel processing, multi-threaded linear algebra libraries will just slow things down. If you have installed OpenBLAS you may need to disable multi-threading. You can see the BLAS R is using in `sessionInfo()`. Disable multi-threading using the `RhpcBLASctl` package:

```
RhpcBLASctl::blas_set_num_threads(1)
```

This needs to be done before using `BiocParallel::bpup`. In the default case of using `MulticoreParam` and not having used `bpup`, weitrix temporarily start a worker pool for large computations, and ensure this is set for workers. If you stray from this default case, we assume you vaguely know what you are doing.




