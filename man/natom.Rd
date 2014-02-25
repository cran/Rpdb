\name{natom}
\alias{natom}
\title{
Number of Atoms in an Object Containing Atomic Coordinates
}
\description{
Evaluates the number of atoms in an object containing atomic coordinates.
}
\usage{
natom(x, ...)
}
\arguments{
  \item{x}{an R object used to select a method.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{Generic function to evalute the number of atom in an object containing atomic coordinates.}
\value{Return an integer indicating the number of atom in \code{x}.}

\seealso{
\code{\link{natom.coords}}, \code{\link{coords}}
}
\examples{
n <- 2
x <- coords(1:n, 1:n, 1:n)
natom(x)
n <- 3
x <- coords(1:n, 1:n, 1:n)
natom(x)

}
\keyword{ array }
