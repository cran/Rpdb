\name{is.conect}
\alias{is.conect}
\title{
Is an Object of Class \sQuote{conect}?
}
\description{
Checks whether its argument is an object of class \sQuote{conect}.
}
\usage{
is.conect(x)
}
\arguments{
  \item{x}{an R object to be tested.}
}
\details{
Tests if x is an object of class \sQuote{conect}, i.e. if x has a \dQuote{class} attribute equal to \code{conect}.
}
\value{
TRUE if x is an object of class \sQuote{conect} and FALSE otherwise
}

\seealso{
\code{\link{conect}}, \code{\link{pdb}}
}
\examples{
x <- conect(eleid.1=rep(1,4), eleid.2=2:5)
is.conect(x)
is.conect(1)

}
\keyword{ classes }
