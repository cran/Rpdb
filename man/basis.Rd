\name{basis}
\alias{basis}
\alias{basis<-}
\alias{basis.default}
\alias{basis<-.default}
\title{
The Basis of an Object
}
\description{
Functions to get or set the basis of an object containing atomic coordinates.
}
\usage{
basis(x)
basis(x) <- value

\method{basis}{default}(x)
\method{basis}{default}(x) <- value

}
\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{value}{a single element character vector use to set the basis of \code{x}.}
}
\details{
\code{basis} and \code{basis<-} are respectively generic accessor and replacement functions.
The default methods get and set the \code{basis} attribute of an object containing atomic coordinates.
This attribute indicate the type basis vector used to express atomic coordinates.
\code{value} must be equal to \code{"xyz"}, for Cartesian, or \code{"abc"}, for fractional coordinates.

}
\value{
\describe{
\item{For \code{basis}:}{NULL or a single element character vector. (NULL is given if the object has no \code{basis} attribute.)}
\item{For \code{basis<-}:}{the updated object.
(Note that the value of \code{basis(x) <- value} is that of the assignment, value, not the return value from the left-hand side.)}
}
}

\seealso{
\code{\link{coords}}, \code{\link{atoms}}
}
\examples{
x <- coords(1:3, 4:6, 7:9, basis= "xyz")
basis(x)
basis(x) <- "abc"
basis(x)
}
\keyword{ attribute }
