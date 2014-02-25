\name{mirrorHelpers}
\alias{Mxy}
\alias{Mxy.coords}
\alias{Mxy.pdb}
\alias{Myz}
\alias{Myz.coords}
\alias{Myz.pdb}
\alias{Mzx}
\alias{Mzx.coords}
\alias{Mzx.pdb}
\alias{Mab}
\alias{Mab.coords}
\alias{Mab.pdb}
\alias{Mbc}
\alias{Mbc.coords}
\alias{Mbc.pdb}
\alias{Mca}
\alias{Mca.coords}
\alias{Mca.pdb}
\title{Helper Functions for reflection of Atomic Coordinates}
\description{Reflection of atomic coordinates with respect to a specific Cartesian plan or a plan defined by two lattice vectors.}
\usage{
Mxy(...)
Myz(...)
Mzx(...)
Mab(...)
Mbc(...)
Mca(...)

\method{Mxy}{coords}(x, mask = TRUE, cryst1 = NULL, \dots)
\method{Myz}{coords}(x, mask = TRUE, cryst1 = NULL, \dots)
\method{Mzx}{coords}(x, mask = TRUE, cryst1 = NULL, \dots)
\method{Mab}{coords}(x, cryst1, mask = TRUE, \dots)
\method{Mbc}{coords}(x, cryst1, mask = TRUE, \dots)
\method{Mca}{coords}(x, cryst1, mask = TRUE, \dots)

\method{Mxy}{pdb}(x, mask = TRUE, cryst1 = x$cryst1, \dots)
\method{Myz}{pdb}(x, mask = TRUE, cryst1 = x$cryst1, \dots)
\method{Mzx}{pdb}(x, mask = TRUE, cryst1 = x$cryst1, \dots)
\method{Mab}{pdb}(x, cryst1 = x$cryst1, mask = TRUE, \dots)
\method{Mbc}{pdb}(x, cryst1 = x$cryst1, mask = TRUE, \dots)
\method{Mca}{pdb}(x, cryst1 = x$cryst1, mask = TRUE, \dots)
}
\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{mask}{a logical vector indicating the set of coordinates to which to apply the reflection.}
  \item{cryst1}{an object of class \sQuote{cryst1} use to convert fractional into Cartesian coordinates when need.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
These functions are helper functions to perform a reflection with respect to a specific Cartesian plan or a plan defined by two lattice vectors. All of them call the \code{mirror} function.
}
\value{An object of the same class as \code{x} with reflected coordinates.}

\seealso{
\code{\link{mirror}} and \code{\link{xyz2abc}}, \code{\link{abc2xyz}} for passing from Cartesian to fractional coordinates (or Vis Versa).

}
\examples{
# First lets read a pdb file
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
visualize(x)
# Mirror operation with respect to the ab-plan
visualize(Mab(x))
# Mirror operation with respect to the ab-plan for residue 1
visualize(Mab(x, mask=x$atoms$resid==1))
}
\keyword{ manip }
