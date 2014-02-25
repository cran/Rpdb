\name{cell.coords}
\alias{cell.coords}
\alias{cell.coords.default}
\alias{cell.coords.cryst1}
\alias{cell.coords.pdb}
\title{
Unit Cell Vectors Coordinates
}
\description{
Computes the Cartesian coordinates of lattice vectors.
}
\usage{
cell.coords(\dots)

\method{cell.coords}{default}(abc, abg = c(90, 90, 90), digits = 3, \dots)

\method{cell.coords}{cryst1}(x, digits = 3, \dots)

\method{cell.coords}{pdb}(x, digits = 3, \dots)
}
\arguments{
  \item{\dots}{arguments passed to methods.}
  \item{abc}{a length 3 numeric vector containing the length of the a, b and c lattice vectors.}
  \item{abg}{a length 3 numeric vector containing the angles (degrees) between the a, b and c lattice vectors (alpha, beta, gamma).}
  \item{digits}{an integer used to round the returned matrix.}
  \item{x}{an R object containing lattice parameters.}
}
\details{
These functions compute a 3x3 matrix whose columns contrain the Cartesian coordinates of lattice vectors.
The 'a' and 'b' vectors are assumed to be respectively along the x-axis and in the xy-plane.

The default method takes directly the lattice parameters as arguments.
For objects of class \code{\link{cryst1}} the lattice parameters are first extracted from the object and then the default method is called.
For objects of class \code{\link{pdb}} the lattice parameters are extracted from their \code{cryst1} component and the default method is called.
}
\value{
A 3x3 matrix containing the Cartesian coordinates of lattice vectors arranged by columns.
}

\seealso{
\code{\link{cryst1}}, \code{\link{xyz2abc}}, \code{\link{pdb}}
}
\examples{
x <- cell.coords(abc = 1:3)
print(x)
x <- cell.coords(abc = 1:3, abg = c(90, 90, 60))
print(x)
c <- cryst1(abc = 1:3, abg = c(90, 90, 60))
is(c)
x <- cell.coords(c)
print(x)
}
\keyword{ manip }
