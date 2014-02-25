\name{coords.pdb}
\alias{coords.pdb}
\alias{coords<-.pdb}
\title{
The Coordinates of an Object of Class \sQuote{pdb}
}
\description{
Get or set the coordinates of an object of class \sQuote{pdb}. 
}
\usage{
\method{coords}{pdb}(x, \dots)

\method{coords}{pdb}(x) <- value
}
\arguments{
  \item{x}{an object of class \sQuote{pdb}.}
  \item{\dots}{further arguments passed to or from other methods.}
  \item{value}{an object of class \sQuote{coords}}
}
\details{
The \code{coords} function uses the \code{atoms} component of a \sQuote{pdb} object to create an object of class \sQuote{coords}.
The \code{coords<-} function sets the atomic coordinates of the \code{atoms} component of a \sQuote{pdb} object. (See \code{\link{coords.atoms}})
}
\value{
The getter function returns a data.frame of class \code{coords} whose columns contain the atomic coordinates of the atoms of a molecular system.
The coordinates can either be Cartesian or fractional coordinates depending on the 'basis' attribute of \code{x}.

The setter function return a \code{pdb} object with updated coordinates.
}


\seealso{
\code{\link{basis}}, \code{\link{coords}}, \code{\link{atoms}}, \code{\link{coords.atoms}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
n <- natom(x)
coords(x) <- coords(rnorm(n), rnorm(n), rnorm(n))
print(x)
class(x)
x <- coords(x)
print(x)
class(x)

}
\keyword{ classes }
\keyword{ manip }
