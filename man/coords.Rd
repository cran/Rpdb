\name{coords}
\alias{coords}
\alias{coords.default}
\alias{coords<-}
\title{
The Atomic Coordinates of an Object
}
\description{
Get or set the atomic coordinates (either Cartesian or fractional coordinates) of an object.
}
\usage{
coords(...)

coords(x) <- value

\method{coords}{default}(x1, x2, x3, basis = "xyz", \dots)
}
\arguments{
  \item{\dots}{arguments passed to methods.}
  \item{x1, x2, x3}{a numeric vector containing the first, second and third coordinates.}
  \item{basis}{a single element character vector indicating the type of basis vector used to express the atomic coordinates.}
  \item{x}{an R object containing atomic coordinates.}
  \item{value}{an object of class \sQuote{coords}.}
}
\details{
The purpose of the \sQuote{coords} class is to store the coordinates of a molecular system
and facilitate their manipulation when passing from the Cartesian to fractional coordinates and vice versa.
\code{coords} and \code{coords<-} are generic accessor and replacement functions.
The default method of the accessor function is actually a builder allowing to create \sQuote{coords} objects from its different components, i.e.:
\code{x1}, \code{x2}, \code{x3}, and \code{basis}.
All the arguments have to be specified except 'basis' which by default is set to "xyz" (Cartesian coordinates).
}
\value{
Return a data.frame of class \sQuote{coords} whose columns contain the three coordinates of the atoms of a molecular system.
The coordinates can either be Cartesian (\code{basis} attribute equal to \code{"xyz"})
or fractional coordinates (\code{basis} attribute equal to \code{"abc"}).
}

\seealso{
\code{\link{basis}}
}
\examples{
x <- coords(c(0,0),c(0,0),c(0,1), basis = "xyz")
print(x)
class(x)
cell <- cryst1(abc = c(10,10,10))
y <- xyz2abc(x,cryst1=cell)
print(y)
class(y)
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
coords(x)
coords(x$atoms)
coords(x) <- coords(x$atoms)

}
\keyword{ classes }
\keyword{ manip }
