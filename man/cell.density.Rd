\name{cell.density}
\alias{cell.density}
\alias{cell.density.default}
\alias{cell.density.pdb}
\title{Density of a Unit Cell}
\description{Compute the density of a unit cell.}
\usage{
cell.density(...)

\method{cell.density}{default}(masses, volume, ...)

\method{cell.density}{pdb}(x, ...)
}
\arguments{
  \item{\dots}{further arguments passed to or from other methods.}
  \item{masses}{a numeric vector containing atomic masses.}
  \item{volume}{a single element numeric vector containing the volume of the unit cell in Angstrom cube.}
  \item{x}{an object of class \sQuote{pdb} containing a \code{cryst1} componant.}
}
\details{
\code{cell.density} is a generic function to compute the density of a unit cell.

For objects of class \sQuote{pdb}:
\itemize{
\item First the volume of the unit cell is calculated by calling the \code{cell.volume} function on the \code{cryst1} component of the \sQuote{pdb} object.
\item Then the element names are converted into element symbols using the \code{toSymbols} function and their masses are taken from the \code{elements} data set.
\item Finally the density is calculated using the sum of the atomic masses and the volume of the unit cell.
}
}
\value{Return a single element numeric vector containing the density of the unit cell in g.cm-3.}

\seealso{
\code{\link{cryst1}}, \code{\link{pdb}}, \code{\link{cell.volume}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
cell.density(x)
}
\keyword{ manip }
