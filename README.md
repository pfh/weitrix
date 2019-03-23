# weitrix

(under development) R package for working with matrices with precision weights.

Install with:

```
BiocManager::install("pfh/weitrix")
```

* [GitHub page](https://github.com/pfh/weitrix)

## Weights

A weight is one over the variance of a measurement. (There may possibly be a further scaling factor. Weights produced by the weitrix package will aim for this to not be needed, so the weight is directly interpretable as a precision.)

A weight of 0 indicates missing data.

## Blessed be the weitrix

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

The `$col` matrix of a `Components` may be used as a design matrix. **Warning:** This may produce liberal results, because the design matrix is itself uncertain and this isn't taken into account. Use this with caution.

