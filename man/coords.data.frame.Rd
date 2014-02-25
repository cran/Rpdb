\name{coords.data.frame}
\alias{coords.data.frame}
\alias{coords.matrix}
\title{
Create \sQuote{coords} Object From \sQuote{data.frame} or \sQuote{matrix}
}
\description{
Create an object of class \sQuote{coords} from a \sQuote{data.frame} or a \sQuote{matrix}.
}
\usage{
\method{coords}{data.frame}(x, basis = NULL, \dots)

\method{coords}{matrix}(x, basis = NULL, \dots)
}
\arguments{
  \item{x}{a data.frame containing x, y and z or a, b and c columns.}
  \item{basis}{a single element character vector indicating the type of basis vector used to express the atomic coordinates.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
When \code{basis==NULL} this function search x, y, z or a, b, c columns in \code{x}.

If x, y, z columns are found in \code{x} they are used to a set the first, second and third coordinates of the returned \sQuote{coords} object.
In that case the basis set of \code{x} is set to \code{"xyz"}.

If a, b, c columns are found in \code{x} they are used to a set the first, second and third coordinates of the returned \sQuote{coords} object.
In that case the basis set of \code{x} is set to \code{"abc"}.

If the function doesn't found neither the x, y, z nor the a, b, c columns an error is returned.

When \code{basis!=NULL} it has to be equal to \code{"xyz"} or \code{"abc"} and \code{x} must have exactly 3 columns.
}
\value{
Return a data.frame of class \sQuote{coords} whose columns contain the three coordinates of the atoms of a molecular system.
The coordinates can either be Cartesian (\code{basis} attribute equal to \code{"xyz"}) or fractional coordinates (\code{basis} attirbute equal to \code{"abc"}).
}

\seealso{
\code{\link{basis}}, \code{\link{coords}}
}

\examples{
x <- data.frame(x=1:3,y=4:6,z=7:9)
class(x)
x <- coords(x)
class(x)
x <- data.frame(a=1:3,b=4:6,c=7:9)
class(x)
x <- coords(x)
class(x)

\dontrun{
x <- data.frame(x=1:3,a=4:6,V=7:9)
x <- coords(x)
}

}
\keyword{ manip }
