\name{is.coords}
\alias{is.coords}
\title{
Is an Object of Class \sQuote{coords}?
}
\description{
Checks whether its argument is an object of class \sQuote{coords}.
}
\usage{
is.coords(x)
}
\arguments{
  \item{x}{
an R object to be tested
}
}
\details{
Tests if x is an object of class \sQuote{coords}, i.e. if x has a \dQuote{class} attribute equal to \code{coords}.
}
\value{
TRUE if x is an object of class \sQuote{coords} and FALSE otherwise
}

\seealso{
\code{\link{coords}}, \code{\link{pdb}}
}
\examples{
x <- coords(c(0,0),c(0,0),c(0,1), basis = "xyz")
is.coords(x)
is.coords(1)

}
\keyword{ classes }
