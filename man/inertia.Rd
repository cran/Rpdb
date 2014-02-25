\name{inertia}
\alias{inertia}
\alias{inertia.coords}
\alias{inertia.atoms}
\alias{inertia.pdb}
\title{Moment of Inertia of a Molecular System}
\description{Computes the inertia tensor of a molecular system from atomic coordinates and masses.}
\usage{
inertia(...)
\method{inertia}{coords}(x, m, \dots)
\method{inertia}{atoms}(x, m = NULL, \dots)
\method{inertia}{pdb}(x, m = NULL, \dots)
}
\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{m}{a numeric vector containing atomic masses.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{\code{inertia} is a generic function to compute the inertia tensor of a molecular system. For object of class \sQuote{coords} both atomic coordinates and masses have to be speifyed. For object of class \sQuote{atoms} the masses are determined from the \code{elename} component of the object (see \code{\link{toSymbols}} and \code{\link{masses}}). For object of class \sQuote{pdb} the \code{atoms} component is used.
}
\value{Return the inertia tensor in a 3x3 matrix.}

\seealso{
\code{\link{toSymbols}}, \code{\link{masses}}, \code{\link{viewInertia}}
}
\examples{
C70 <- read.pdb(system.file("examples/C70.pdb",package="Rpdb"))
inertia(C70)
visualize(C70, mode = NULL)
viewXY()
viewInertia(C70)
}
\keyword{ manip }
