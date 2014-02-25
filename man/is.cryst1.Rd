\name{is.cryst1}
\alias{is.cryst1}
\title{
Is an Object of Class \sQuote{cryst1}?
}
\description{
Checks whether its argument is an object of class \sQuote{cryst1}.
}
\usage{
is.cryst1(x)
}
\arguments{
  \item{x}{
an R object to be tested
}
}
\details{
Tests if x is an object of class \sQuote{cryst1}, i.e. if x has a \dQuote{class} attribute equal to \code{cryst1}.
}
\value{
TRUE if x is an object of class \sQuote{cryst1} and FALSE otherwise
}

\seealso{
\code{\link{cryst1}}, \code{\link{pdb}}
}
\examples{
x <- cryst1(abc=rep(10,3))
is.cryst1(x)
is.cryst1(1)

}
\keyword{ classes }
