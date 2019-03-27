
#' @import methods
#' @import utils
#' @import stats
#'
#' @importFrom assertthat
#'     assert_that is.string is.number
#'
#' @importFrom purrr
#'     map map_chr map_dbl map_int map2
#'
#' @importFrom limma
#'     getEAWP
#'
#' @importClassesFrom SummarizedExperiment 
#'     SummarizedExperiment
#' @importFrom SummarizedExperiment
#'     SummarizedExperiment assay assay<- rowData colData
#'
#' @importFrom DelayedArray
#'     realize getAutoBPPARAM getAutoBlockSize
#'
#' @importFrom DelayedMatrixStats
#'     rowSums2 colSums2
#'
#' @importFrom S4Vectors
#'     metadata metadata<-
#'
#' @importFrom BiocParallel
#'     bplapply
#'
#' @importFrom BiocGenerics
#'     rowSums colSums rbind cbind
#'
NULL

