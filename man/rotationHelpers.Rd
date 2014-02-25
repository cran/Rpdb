\name{rotationHelpers}
\alias{Rx}
\alias{Rx.coords}
\alias{Rx.pdb}
\alias{Ry}
\alias{Ry.coords}
\alias{Ry.pdb}
\alias{Rz}
\alias{Rz.coords}
\alias{Rz.pdb}
\title{Helper Functions for Rotation of Atomic Coordinates}
\description{Rotation of atomic coordinates along a specific Cartesian vector.}
\usage{
Rx(...)
Ry(...)
Rz(...)

\method{Rx}{coords}(x, angle = 0, mask = TRUE, cryst1 = NULL, \dots)
\method{Ry}{coords}(x, angle = 0, mask = TRUE, cryst1 = NULL, \dots)
\method{Rz}{coords}(x, angle = 0, mask = TRUE, cryst1 = NULL, \dots)

\method{Rx}{pdb}(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, \dots)
\method{Ry}{pdb}(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, \dots)
\method{Rz}{pdb}(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, \dots)
}
\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{angle}{the angle of the rotation.}
  \item{mask}{a logical vector indicating the set of coordinates to which the rotation has to be applyed.}
  \item{cryst1}{an object of class \sQuote{cryst1} use to convert fractional into Cartesian coordinates when need.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
These functions are helper functions to perform a rotation around a specific Cartesian vector. All of them call the \code{R} function.
}
\value{An object of the same class as \code{x} with rotated coordinates.}
\seealso{
\code{\link{R}} and \code{\link{xyz2abc}}, \code{\link{abc2xyz}} for passing from Cartesian to fractional coordinates (or Vis Versa).
}
\examples{
# First lets read a pdb file
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
cell <- cell.coords(x)
visualize(x)
# Rotation of the structure around the z-axis
visualize(Rz(x, 90))
# Rotation of the residue 1 around the c-axis
visualize(Rz(x, 90, mask=x$atoms$resid==1))
}
\keyword{ manip }
