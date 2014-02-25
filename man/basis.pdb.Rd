\name{basis.pdb}
\alias{basis.pdb}
\alias{basis<-.pdb}
\title{
The Basis Set of a \sQuote{pdb} Object}
\description{
Functions to get or set the basis of a 'pdb' object.
}
\usage{
\method{basis}{pdb}(x)

\method{basis}{pdb}(x) <- value
}
\arguments{
  \item{x}{an object of class \sQuote{pdb}.}
  \item{value}{a single element character vector used to set the basis of \code{x}.}
}
\details{
These two function get and set the \code{\link{basis}} attribute of the \code{\link{atoms}} component of the \code{\link{pdb}} object.
This attribute indicate the type basis vector used to express atomic coordinates.
\code{value} must be equal to \code{"xyz"}, for Cartesian, or \code{"abc"}, for fractional coordinates.

}
\value{
\describe{
\item{For \code{basis}:}{a single element character vector equal to \code{"xyz"} or \code{"abc"} for Cartesian or fractional coordinates respectively.}
\item{For \code{basis<-}:}{the updated \sQuote{pdb} object with the new basis.}
}
}

\seealso{
\code{\link{basis}}, \code{\link{pdb}}
}

\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
basis(x)
basis(x) <- "abc"
basis(x)

}
\keyword{ attribute }