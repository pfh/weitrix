
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
#'     map map_chr map_dbl map_int map_lgl map2
#'
#' @importFrom ggplot2
#'    ggplot aes aes_string geom_point geom_line geom_hline labs coord_cartesian
#'    geom_smooth coord_trans geom_hline geom_boxplot geom_abline theme
#'    element_text element_blank facet_wrap vars geom_bin2d scale_fill_viridis_c
#'
#' @importFrom rlang
#'    eval_tidy enquo quo quo_get_expr as_label .data
#'
#' @importFrom scales
#'    percent
#'
#' @importFrom limma
#'     getEAWP squeezeVar
#'
#' @importFrom topconfects
#'     nest_confects normal_confects
#'
#' @importClassesFrom SummarizedExperiment 
#'     SummarizedExperiment
#' @importFrom SummarizedExperiment
#'     SummarizedExperiment assay assay<- rowData rowData<- colData
#'
#' @importFrom DelayedArray
#'     realize getAutoBPPARAM setAutoBPPARAM getAutoBlockSize 
#'     sweep apply
#'
#' @importFrom DelayedMatrixStats
#'     rowSums2 colSums2
#'
#' @importFrom S4Vectors
#'     metadata metadata<-
#'
#' @importFrom BiocParallel
#'     bplapply bpnworkers bpisup bpnworkers bpstart bpstop bpparam
#'     SerialParam
#'
#' @importFrom RhpcBLASctl
#'     blas_get_num_procs blas_set_num_threads
#'
#
# Missing from Bioconductor 3.17 2023-04-07. Let's hope it works without it.
# @importFrom BiocGenerics
#     rowSums colSums
#
#' @importFrom BiocGenerics
#'     rbind cbind
#'
#' @importFrom splines
#'     ns
#'
#' @importFrom Ckmeans.1d.dp
#'     Ckmeans.1d.dp
#'
#' @importFrom glm2
#'     glm2
#'
#' @importFrom grDevices
#'     dev.cur
#'
NULL

