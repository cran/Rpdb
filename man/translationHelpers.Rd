\name{translationHelpers}
\alias{Tx}
\alias{Tx.coords}
\alias{Tx.pdb}
\alias{Ty}
\alias{Ty.coords}
\alias{Ty.pdb}
\alias{Tz}
\alias{Tz.coords}
\alias{Tz.pdb}
\alias{Ta}
\alias{Ta.coords}
\alias{Ta.pdb}
\alias{Tb}
\alias{Tb.coords}
\alias{Tb.pdb}
\alias{Tc}
\alias{Tc.coords}
\alias{Tc.pdb}
\title{Helper Functions for Translation of Atomic Coordinates}
\description{Translation of atomic coordinates along a specific Cartesian or lattice vector.}
\usage{
### Tranlastion of Cartesian coordinates:
Tx(...)
Ty(...)
Tz(...)

\method{Tx}{coords}(obj, x = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, \dots)
\method{Ty}{coords}(obj, y = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, \dots)
\method{Tz}{coords}(obj, z = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, \dots)

\method{Tx}{pdb}(obj, x = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, \dots)
\method{Ty}{pdb}(obj, y = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, \dots)
\method{Tz}{pdb}(obj, z = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, \dots)

### Tranlastion of fractional coordinates:
Ta(...)
Tb(...)
Tc(...)

\method{Ta}{coords}(obj, a = 0, mask = TRUE, cryst1 = NULL, \dots)
\method{Tb}{coords}(obj, b = 0, mask = TRUE, cryst1 = NULL, \dots)
\method{Tc}{coords}(obj, c = 0, mask = TRUE, cryst1 = NULL, \dots)

\method{Ta}{pdb}(obj, a = 0, mask = TRUE, cryst1 = obj$cryst1, \dots)
\method{Tb}{pdb}(obj, b = 0, mask = TRUE, cryst1 = obj$cryst1, \dots)
\method{Tc}{pdb}(obj, c = 0, mask = TRUE, cryst1 = obj$cryst1, \dots)
}
\arguments{
  \item{obj}{an R object containing atomic coordinates.}
  \item{x}{the x-component of the translation vector.}
  \item{y}{the y-component of the translation vector.}
  \item{z}{the z-component of the translation vector.}
  \item{a}{the a-component of the translation vector.}
  \item{b}{the b-component of the translation vector.}
  \item{c}{the c-component of the translation vector.}
  \item{mask}{a logical vector indicating the set of coordinates to which to apply the translation.}
  \item{thickness}{a numeric value indicating the fraction of the thicknees of the selected atom to be added to the translation vector (Usually 0, 0.5 or 1. See details).}
  \item{cryst1}{an object of class \sQuote{cryst1} use to convert Cartesian into fraction coordinates (or Vis Versa) when need.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
These functions are helper functions to perform a translation along a specific Cartesian or lattice vector. All of them call either the \code{Txyz} or \code{Tabc} function.
}
\value{An object of the same class as \code{x} with translated coordinates.}

\seealso{
\code{\link{Txyz}}, \code{\link{Tabc}}\cr
Passing from Cartesian to fractional coordinates (or Vis Versa):\cr
\code{\link{xyz2abc}}, \code{\link{abc2xyz}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
visualize(x, mode = NULL)
visualize(Ty(x, 10), mode = NULL)
visualize(Ty(x, 10, mask=x$atoms$resid==1), mode = NULL)
visualize(Tb(x, 1 ), mode = NULL)
visualize(Tb(x, 1 , mask=x$atoms$resid==1), mode = NULL)

# Lets build a C70/Pentacene dimer with an inter-molecular distance equal to 3.5
C70 <- read.pdb(system.file("examples/C70.pdb",package="Rpdb"))
Pen <- read.pdb(system.file("examples/Pentacene.pdb",package="Rpdb"))
x <- merge(C70,Pen)
visualize(x, mode = NULL)
viewXY()
visualize(Tz(x, z=3.5, mask=x$atoms$resname=="C70", thickness=0.5), mode = NULL)
viewXY()
}
\keyword{ manip }
