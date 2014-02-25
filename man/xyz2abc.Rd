\name{xyz2abc}
\alias{abc2xyz}
\alias{xyz2abc}
\alias{abc2xyz.coords}
\alias{xyz2abc.coords}
\alias{abc2xyz.atoms}
\alias{xyz2abc.atoms}
\alias{abc2xyz.pdb}
\alias{xyz2abc.pdb}
\title{
From Cartesian to Fractional Coordinates and Vis Versa
}
\description{
Converts Cartesian coordinates into fractional coordinates and vice versa.
}
\usage{
xyz2abc(\dots)
abc2xyz(\dots)

\method{xyz2abc}{coords}(x, cryst1, \dots)
\method{abc2xyz}{coords}(x, cryst1, \dots)

\method{xyz2abc}{atoms}(x, cryst1, \dots)
\method{abc2xyz}{atoms}(x, cryst1, \dots)

\method{xyz2abc}{pdb}(x, cryst1 = x$cryst1, \dots)
\method{abc2xyz}{pdb}(x, cryst1 = x$cryst1, \dots)
}
\arguments{
  \item{\dots}{arguments passed to methods.}
  \item{x}{an R object containing atomic coordinates.}
  \item{cryst1}{an object of class \code{\link{cryst1}}.}
}
\details{
For \code{\link{atoms}} and \code{\link{pdb}} objects, the atomic coordinates
are first extracted from \code{x} using the \code{\link{coords}} function.
Then, using the periodic boundary conditions stored into \code{cryst1}, the
coordinates are converted from Cartesian to fractional (for the \code{xyz2abc} functions)
or from fractional to Cartesian (for the \code{abc2xyz} functions) coordinates.
Finally, for \code{\link{atoms}} and \code{\link{pdb}} objects, the new atomic coordinates
are reassigned to the original \code{x} object using the \code{\link{coords<-}} function and \code{x} is returned.
}
\value{
Return an object of the same class as \code{x}, with atomic coordinates expressed in a different basis set.
}

\seealso{
\code{\link{basis}}, \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}, \code{\link{cryst1}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
basis(x)
x <- xyz2abc(x)
basis(x)
x <- abc2xyz(x)
basis(x)

\dontrun{

# This example return an error because the coordinates stored
# into the PDB file are already Cartesian coordinates.
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
x <- abc2xyz(x)
}
}
\keyword{ manip }
