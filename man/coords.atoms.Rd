\name{coords.atoms}
\alias{coords.atoms}
\alias{coords<-.atoms}
\title{
The Atomic Coordinates of an Object of Class \sQuote{atoms}
}
\description{
Get or set the atomic coordinates of an object of class \sQuote{atoms}.
}
\usage{
\method{coords}{atoms}(x, \dots)

\method{coords}{atoms}(x) <- value
}
\arguments{
  \item{x}{an object of class \sQuote{atoms}.}
  \item{\dots}{further arguments passed to or from other methods.}
  \item{value}{an object of class \sQuote{coords}.}
}
\details{
The \code{coords} function extract the \code{x1}, \code{x2} and \code{x3} components of a \sQuote{atoms} object
as well as its \code{basis} attribute to create a \sQuote{coords} object.
The \code{coords<-} function set the \code{x1}, \code{x2} and \code{x3} components of a \sQuote{atoms} object as well as its \code{basis} attribute.
}
\value{
The getter function return a data.frame of class \code{coords} whose columns contain the atomic coordinates of each element of a molecular system.
The coordinates can either be Cartesian or fractional coordinates depending on the \code{basis} attribute of \code{x}.

The setter function return a \code{atoms} object with updated atomic coordinates.
}

\seealso{
\code{\link{basis}}, \code{\link{coords}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
x <- x$atoms
n <- natom(x)
coords(x) <- coords(rnorm(n), rnorm(n), rnorm(n))
head(x)
class(x)
x <- coords(x)
head(x)
class(x)

}
\keyword{ classes }
\keyword{ manip }
