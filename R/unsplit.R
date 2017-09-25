#' Reassemble Groups
#' 
#' \code{unsplit} reverses the effect of \code{split}.
#' 
#' \code{unsplit} is a generic functions with a default method (Method dispatch 
#' takes place based on the class of the first element of \code{value}) working 
#' with lists of vectors or data frames (assumed to have compatible structure, 
#' as if created by \code{split}). It puts elements or rows back in the 
#' positions given by \code{f}. In the data frame case, row names are obtained 
#' by unsplitting the row name vectors from the elements of \code{value}. \cr 
#' \code{f} is recycled as necessary and if the length of \code{x} is not a 
#' multiple of the length of \code{f} a warning is printed. \cr Any missing 
#' values in \code{f} are dropped together with the corresponding values of 
#' \code{x}.
#' 
#' @return Returns a vector or data frame for which \code{split(x, f)} equals
#'   \code{value}
#'   
#' @param value a list of vectors or data frames compatible with a splitting of
#'   \code{x}. Recycling applies if the lengths do not match.
#' @param f a \sQuote{factor} in the sense that \code{\link{as.factor}(f)}
#'   defines the grouping, or a list of such factors in which case their
#'   interaction is used for the grouping.
#' @param drop logical indicating if levels that do not occur should be dropped
#'   (if \code{f} is a \code{factor} or a list).
#' @param \dots further potential arguments passed to methods.
#' 
#' @seealso 
#' \code{\link[base]{cut}} to categorize numeric values.
#' \cr
#' \code{\link[base]{strsplit}} to split strings.
#' 
#' @references
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
#' 
#' @examples
#' require(stats); require(graphics)
#' n <- 10; nn <- 100
#' g <- factor(round(n * runif(n * nn)))
#' x <- rnorm(n * nn) + sqrt(as.numeric(g))
#' xg <- split(x, g)
#' boxplot(xg, col = "lavender", notch = TRUE, varwidth = TRUE)
#' sapply(xg, length)
#' sapply(xg, mean)
#' 
#' ### Calculate 'z-scores' by group (standardize to mean zero, variance one)
#' z <- unsplit(lapply(split(x, g), scale), g)
#' 
#' # or
#' 
#' zz <- x
#' split(zz, g) <- lapply(split(x, g), scale)
#' 
#' # and check that the within-group std dev is indeed one
#' tapply(z, g, sd)
#' tapply(zz, g, sd)
#' 
#' 
#' ### data frame variation
#' 
#' ## Notice that assignment form is not used since a variable is being added
#' 
#' g <- airquality$Month
#' l <- split(airquality, g)
#' l <- lapply(l, transform, Oz.Z = scale(Ozone))
#' aq2 <- unsplit(l, g)
#' head(aq2)
#' with(aq2, tapply(Oz.Z,  Month, sd, na.rm=TRUE))
#' 
#' 
#' ### Split a matrix into a list by columns
#' ma <- cbind(x = 1:10, y = (-4:5)^2)
#' split(ma, col(ma))
#' 
#' split(1:10, 1:2)
#'
#' @keywords category
#' 
#' @name unsplit
#' @export        
unsplit <- function(value, f, drop = FALSE, ...)
  UseMethod("unsplit", value[[1]])

#' @rdname unsplit
#' @export
unsplit.default <- function(value, f, drop = FALSE, ...)
{
  len <- length(if (is.list(f)) f[[1L]] else f)
  if (is.data.frame(value[[1L]])) {
    x <- value[[1L]][rep(NA, len), , drop = FALSE]
    rownames(x) <- unsplit(lapply(value, rownames), f, drop = drop)
  }
  else x <- value[[1L]][rep(NA, len)]
  split(x, f, drop = drop) <- value
  x
}
