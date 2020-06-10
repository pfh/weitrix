
#' Convert data to a weitrix
#'
#' Ensure data is a weighted matrix or "weitrix". A weitrix is a
#' SummarizedExperiment or subclass thereof with some metadata fields set. If it
#' is ambiguous how to do this, produce an error.
#'
#' Input can be a matrix or DelayedArray.
#'
#' Input can be anything the limma package recognizes, notably the limma EList
#' class (for example as output by \code{voom} or \code{vooma}).
#'
#' If weights are not present in "object" and not given with "weights", they
#' default for 0 for NA values and 1 for everything else.
#'
#' @param object Object to convert.
#' @param weights Optional, weights matrix if not present in \code{object}.
#'
#' @return A SummarizedExperiment object with metadata fields marking it as a
#' weitrix.
#'
#' @examples
#' mat <- matrix(c(1,2,NA,3,NA,4), ncol=2)
#' weitrix <- as_weitrix(mat)
#'
#' metadata(weitrix)
#' weitrix_x(weitrix)
#' weitrix_weights(weitrix)
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
        if (is.null(weights)) {
            weights <- !is.na(object)
            #Fails: mode(weights) <- "numeric"
            weights <- weights+0
        }

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
        eawp$weights <- as.matrix(weights)
    }

    x <- eawp$exprs
    weights <- eawp$weights
    probes <- eawp$probes

    # Default to weights representing non-missingness
    if (is.null(weights)) {
        weights <- !is.na(x)
        mode(weights) <- "numeric"
    }

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
#' Set metadata entries in a SummarizedExperiment object
#' so that it can be used as a weitrix.
#'
#' @param object A SummarizedExperiment object.
#' @param x_name Name of the assay containing the observations.
#' @param weights_name Name of the assay containing the weights.
#'
#' @return
#' A SummarizedExperiment object with metadata fields marking it as a weitrix.
#'
#' @examples
#' mat <- matrix(c(1,2,NA,3,NA,4), ncol=2)
#' weights <- matrix(c(1,0.5,0,2,0,1), ncol=2)
#' se <- SummarizedExperiment(assays=list(foo=mat, bar=weights))
#' 
#' weitrix <- bless_weitrix(se, "foo", "bar")
#'
#' metadata(weitrix)
#' weitrix_x(weitrix)
#' weitrix_weights(weitrix)
#'
#' @return
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

#' Get or set a weitrix object's "x" matrix 
#'
#' Gets or sets the appropriate assay in the SummarizedExperiment object.
#'
#' @param weitrix A weitrix object.
#'
#' @return 
#' A matrix-like object such as a matrix or a DelayedArray.
#' 
#' @examples
#' weitrix_x(simwei)
#'
#' simwei2 <- simwei
#' weitrix_x(simwei2) <- weitrix_x(simwei2) * 2
#'
#' @export
weitrix_x <- function(weitrix) {
    assay(weitrix, metadata(weitrix)$weitrix$x_name)
}

#' Get or set a weitrix object's "weights" matrix
#'
#' Gets or sets the appropriate assay in the SummarizedExperiment object.
#'
#' @param weitrix A weitrix object.
#'
#' @return 
#' A matrix-like object such as a matrix or a DelayedArray.
#' 
#' @examples
#' weitrix_weights(simwei)
#'
#' @export
weitrix_weights <- function(weitrix) {
    assay(weitrix, metadata(weitrix)$weitrix$weights_name)
}


#' @rdname weitrix_x
#' @param x The weitrix to modify.
#' @param value The new matrix.
#' @export
`weitrix_x<-` <- function(x, value) {
    assay(x, metadata(x)$weitrix$x_name) <- value
    x
}

#' @rdname weitrix_weights
#' @param x The weitrix to modify.
#' @param value The new matrix.
#' @export
`weitrix_weights<-` <- function(x, value) {
    assay(x, metadata(x)$weitrix$weights_name) <- value
    x
}

#' Convert a weitrix object to a limma EList object
#'
#' The resulting object can be used as input to \code{limma::lmFit}
#' for a limma analysis.
#'
#' @param weitrix A weitrix object.
#'
#' @return
#' A limma EList object.
#'
#' @examples
#' library(limma)
#'
#' elist <- weitrix_elist(simwei)
#' design <- model.matrix(~true_score, data=colData(simwei))
#' fit <- lmFit(elist, design)
#' # ...perform further limma analysis steps as desired...
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








