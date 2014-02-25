\name{cryst1}
\alias{cryst1}
\alias{cryst1.default}
\title{
Create \sQuote{cryst1} Object
}
\description{
Create an object of class \sQuote{cryst1} containing the unit cell parameters and the name of the space group to associate with an object of class \sQuote{pdb}.
}
\usage{
cryst1(...)

\method{cryst1}{default}(abc, abg = c(90, 90, 90), sgroup = "P1", \dots)
}
\arguments{
  \item{\dots}{further arguments passed to or from other methods.}
  \item{abc}{a numeric vector of length 3 containing the norms of the lattice vectors a, b and c.}
  \item{abg}{a numeric vector of length 3 containing the angles between the lattice vectors \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.}
  \item{sgroup}{a character string giving the Hermann-Mauguin symbol of the space group.}
}
\details{
This function is the generic function to create objects of class \sQuote{cryst1}.
The purpose of this class is to store CRYST1 records from PDB files which contain the unit cell parameters
and the name of the space group of a molecular system stored in a PDB file.
The default method of the \code{cryst1} function create an object of class \sQuote{cryst1} from its different components,
i.e.: \code{abc}, \code{abg} and \code{sgroup}. At least \code{abc} has to be specified.
}
\value{
Retur a list of class \sQuote{cryst1} with the following components:
\item{abc}{a numeric vector of length 3 containing the norms of the lattice vectors a, b and c.}
\item{abg}{a numeric vector of length 3 containing the angles between the lattice vectors \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.}
\item{sgroup}{a character string giving the Hermann-Mauguin symbol of the space group.}
}

\seealso{
\code{\link{cell.coords}}, \code{\link{pdb}}
}
\examples{
x <- cryst1(abc = c(10, 10, 10), abg = c(90,90,90), sgroup = "P1")
class(x)

}
\keyword{ classes }
