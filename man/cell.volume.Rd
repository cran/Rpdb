\name{cell.volume}
\alias{cell.volume}
\alias{cell.volume.cryst1}
\alias{cell.volume.pdb}
\title{Volume of a Unit Cell}
\description{Computes the volume of a unit cell.}
\usage{
cell.volume(...)

\method{cell.volume}{cryst1}(x, ...)

\method{cell.volume}{pdb}(x, ...)
}
\arguments{
  \item{\dots}{further arguments passed to or from other methods.}
  \item{x}{an R object containing unit cell parameters}
}
\details{
\code{cell.volume} is a generic function to compute the volume of a unit cell.

For objects of class \sQuote{cryst1}, the \code{cell.coords} function is called to build a matrix containing the cartesian coordinates of the lattice vectors. Then, these vector are used to compute the volume of the unit cell.

For objects of class \sQuote{pdb}, their \code{cryst1} component is used.
}
\value{Return a single element numeric vector containing the volume of the unit cell in Angstrom cube.}

\seealso{
\code{\link{cell.coords}}, \code{\link{cryst1}}, \code{\link{pdb}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
cell.volume(x)

}
\keyword{ manip }
