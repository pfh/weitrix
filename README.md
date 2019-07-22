# weitrix

(under development) R package for working with matrices with precision weights.

Install with:

```
BiocManager::install("pfh/weitrix")
```

* [GitHub page](https://github.com/pfh/weitrix)
* [Documentation site](http://logarithmic.net/weitrix/)

## Weights

A weight is one over the variance of a measurement. (There may possibly be a further scaling factor. Weights produced by the weitrix package will aim for this to not be needed, so the weight is directly interpretable as a precision.)

A weight of 0 indicates missing data.

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

weitrix uses DelayedArray assays. weitrix will attempt to perform calculations blockwise in parallel. weitrix tries to use DelayedArray defaults. Adjust with `DelayedArray::setRealizationBackend`, 
`DelayedArray::setAutoBlockSize`, `DelayedArray::setAutoBPPARAM`.

The x and weights assays in a weitrix may be stored in an exotic matrix class, but it should always be possible to convert back to a normal R matrix with `as.matrix`.

Set the DelayedArray realization backend to `RleMatrix` or `HDF5Array` if weitrices will be too big to fit in memory uncompressed. The `RleMatrix` stores data in memory in a compressed form. The `HDF5Array` backend stores data on disk, in temporary files.

* `DelayedArray::setRealizationBackend("RleArray")` may require copying the complete matrix to all members of the cluster (not sure on this point). HDF5Array seems better, since all the actual data is stored on disk and the disk cache will be shared between processes.

* If using `DelayedArray::setRealizationBackend("HDF5Array")` you may also want to set `HDF5Array::setHDF5DumpDir`.

* A weitrix can be permanently stored to disk using `HDF5Array::saveHDF5SummarizedExperiment`.

Example setup:

```
library(DelayedArray)
library(HDF5Array)
library(BiocParallel)

# Start a 4 worker cluster (Linux or Mac OS)
param <- bpstart(MulticoreParam(4))
setAutoBPPARAM(param)

# Store intermediate results in a directory called __dump__
# You may need to clean up this directory manually
setRealizationBackend("HDF5Array")
setHDF5DumpDir("__dump__")
```

## BiocParallel problems

Parallelism in R and Bioconductor remains finicky but is also necessary for large datasets. weitrix uses parallel processing by default.

If weitrix hangs, try running it with serial processing:

```
DelayedArray::setAutoBPPARAM( BiocParallel::SerialParam() )
```




