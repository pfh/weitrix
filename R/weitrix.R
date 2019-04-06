
#' Ensure data is a weighted matrix
#'
#' Ensure data is a weighted matrix or "weitrix". 
#' A weitrix is a SummarizedExperiment or subclass thereof with some metadata fields set. 
#' If it is ambiguous how to do this, produce an error.
#'
#' Input can be a matrix or DelayedArray.
#'
#' Input can be anything the limma package recognizing,
#' notably the EList class (for example as output by voom or vooma).
#'
#' If weights are not present in "object" and not given with "weights",
#' they default for 0 for NA values and 1 for everything else.
#'
#' @export
as_weitrix <- function(object, weights=NULL) {
    if (is(object, "SummarizedExperiment")) {
        assert_that(is.null(weights), msg=
            "Add weights as an assay to the object, then use bless_weights().")
        assert_that(!is.null(metadata(object)$weitrix), msg=paste0(
            "Use bless_weitrix() to specify x and weights assays ",
            "to use in a SummarizedExperiment object."))
        return(object)
    }

    if (is(object, "DelayedArray")) {
        if (is.null(weights))
            weights <- !is.na(object)

        result <- SummarizedExperiment(
            assays=list(x=object, weights=weights))
        result <- bless_weitrix(result, "x", "weights")
        return(result)
    }

    assert_that(requireNamespace("limma", quietly=TRUE), msg=
        "Will try limma::getEAWP to understand object. Please install limma.")

    eawp <- getEAWP(object)
    if (!is.null(weights)) {
        if (!is.null(eawp$weights))
            warning("Overwriting existing weights")
        eawp$weights <- weights
    }

    x <- eawp$exprs
    weights <- eawp$weights
    probes <- eawp$probes

    # Default to weights representing non-missingness
    if (is.null(weights))
        weights <- !is.na(x)

    ## Missingness always encoded as 0 weight, value irrelevant
    assert_that(identical( dim(x), dim(weights) ))
    #assert_that(!any(is.na(weights)))
    #assert_that(!any(is.na(x[ weights != 0 ])))

    result <- SummarizedExperiment(
        assays = list(x = x, weights = weights),
        rowData = probes
    )

    bless_weitrix(result, "x", "weights")
}

#' Bless a SummarizedExperiment as a weitrix
#'
#' @export
bless_weitrix <- function(object, x_name, weights_name) {
    assert_that(is(object, "SummarizedExperiment"))
    assert_that(is.string(x_name))
    assert_that(is.string(weights_name))

    metadata(object)$weitrix <- list(
        x_name = x_name,
        weights_name = weights_name)

    object
}

#' Get a weitrix object's "x" matrix 
#'
#' @export
weitrix_x <- function(weitrix) {
    assay(weitrix, metadata(weitrix)$weitrix$x_name)
}

#' Get a weitrix object's "weights" matrix
#'
#' @export
weitrix_weights <- function(weitrix) {
    assay(weitrix, metadata(weitrix)$weitrix$weights_name)
}

#' @export
`weitrix_x<-` <- function(x, value) {
    assay(x, metadata(x)$weitrix$x_name) <- value
    x
}

#' @export
`weitrix_weights<-` <- function(x, value) {
    assay(x, metadata(x)$weitrix$weights_name) <- value
    x
}

#' Convert a weitrix object to a limma EList object
#'
#' @export
weitrix_elist <- function(weitrix) {
    weitrix <- as_weitrix(weitrix)

    new("EList", list(
        E=as.matrix(weitrix_x(weitrix)),
        weights=as.matrix(weitrix_weights(weitrix)),
        genes=as.data.frame(rowData(weitrix)),
        targets=as.data.frame(colData(weitrix))))
}








