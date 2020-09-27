
#' Natural cubic spline with good choice of knots
#'
#' For use in model formulas, 
#' natural cubic spline as in \code{splines::ns} 
#' but with knot positions chosen using
#' k-means rather than quantiles.
#' Automatically uses less knots if there are insufficient distinct values.
#'
#' Wong (1982, 1984) showed the asymptotic density of k-means in 1 dimension is 
#'     proportional to the cube root of the density of x. 
#' Compared to using quantiles (the default for \code{\link[splines]{ns}}),
#'     choosing knots using k-means produces a better spread of knot locations
#'     if the distribution of values is very uneven.
#'
#' k-means is computed in an optimal, deterministic way using
#' \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}.
#'
#' @param x The predictor variable. A numeric vector.
#'
#' @param n_knots Number of knots to use.
#'
#' @param verbose If TRUE, produce a message about the knots chosen.
#' 
#' @return
#' A matrix of predictors, similar to \code{\link[splines]{ns}}.
#'
#' This function supports "safe prediction" 
#' (see \code{\link[stats]{makepredictcall}}). 
#' Original knot locations will be used for prediction with 
#' \code{\link[stats]{predict}}.
#'
#' @seealso
#' \code{\link[splines]{ns}}, \code{\link[stats]{makepredictcall}}
#'
#' @examples
#' lm(mpg ~ well_knotted_spline(wt,3), data=mtcars)
#'
#' # When insufficient unique values exist, less knots are used
#' lm(mpg ~ well_knotted_spline(gear,3), data=mtcars)
#'
#' library(ggplot2)
#' ggplot(diamonds, aes(carat, price)) + 
#'    geom_point() + 
#'    geom_smooth(method="lm", formula=y~well_knotted_spline(x,10))
#'
#' @references
#'
#' Wong, M. (1982).
#' Asymptotic properties of univariate sample k-means clusters.
#' Working paper #1341-82, Sloan School of Management, MIT.
#' \url{https://dspace.mit.edu/handle/1721.1/46876}
#'
#' Wong, M. (1984). 
#' Asymptotic properties of univariate sample k-means clusters.
#' \emph{Journal of Classification}, 1(1), 255–270.
#' \url{https://doi.org/10.1007/BF01890126}
#'
#' @export
well_knotted_spline <- function(x, n_knots, verbose=TRUE) {
    label <- deparse(substitute(x))
    n_unique <- length(unique(x))
    boundary <- range(x)
    n_knots <- min(n_knots, n_unique-2)
    if (n_knots <= 0) {
        knots <- numeric(0)
    } else {
        knots <- Ckmeans.1d.dp(x, n_knots)$centers
        knots <- knots[! knots %in% boundary ]
    }

    if (verbose)
        message(label, " range ", 
            paste(format(boundary),collapse=" "), 
            " knots ", paste(format(knots),collapse=" "))
    
    spline_function <- quote(splines::ns)
    result <- eval(spline_function)(x, knots=knots, Boundary.knots=boundary)
    class(result) <- c("well_knotted_spline", class(result))
    attributes(result)$spline_function <- spline_function
    result
}


#' @export
makepredictcall.well_knotted_spline <- function(var, call) {
    assert_that(is.call(call))
    call[[1]] <- attributes(var)$spline_function
    class(var) <- class(var)[-1]
    result <- makepredictcall(var, call)
    #print(result)
    result
}
