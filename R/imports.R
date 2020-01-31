
#' @import methods
#' @import utils
#' @import stats
#'
#' @importFrom reshape2
#'     melt
#'
#' @importFrom assertthat
#'     assert_that is.string is.number
#'
#' @importFrom dplyr
#'     left_join
#'
#' @importFrom purrr
#'     map map_chr map_dbl map_int map2
#'
#' @importFrom ggplot2
#'    ggplot aes_string geom_point geom_line geom_hline labs coord_cartesian
#'    
#' @importFrom scales
#'    percent
#'
#' @importFrom limma
#'     getEAWP
#'
#' @importClassesFrom SummarizedExperiment 
#'     SummarizedExperiment
#' @importFrom SummarizedExperiment
#'     SummarizedExperiment assay assay<- rowData rowData<- colData
#'
#' @importFrom DelayedArray
#'     realize getAutoBPPARAM getAutoBlockSize sweep apply
#'
#' @importFrom DelayedMatrixStats
#'     rowSums2 colSums2
#'
#' @importFrom S4Vectors
#'     metadata metadata<-
#'
#' @importFrom BiocParallel
#'     bplapply bpnworkers bpisup bpnworkers bpstart bpstop
#'
#' @importFrom RhpcBLASctl
#'     blas_get_num_procs blas_set_num_threads
#'
#' @importFrom BiocGenerics
#'     rowSums colSums rbind cbind
#'
NULL

