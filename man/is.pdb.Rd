\name{is.pdb}
\alias{is.pdb}
\title{
Is an Object of Class \sQuote{pdb}?
}
\description{
Checks whether its argument is an object of class \sQuote{pdb}.
}
\usage{
is.pdb(x)
}
\arguments{
  \item{x}{
an R object to be tested
}
}
\details{
Tests if x is an object of class \sQuote{pdb}, i.e. if x has a \dQuote{class} attribute equal to \code{pdb}.
}
\value{
TRUE if x is an object of class \sQuote{pdb} and FALSE otherwise
}

\seealso{
\code{\link{pdb}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
is.pdb(x)
is.pdb(1)

}
\keyword{ classes }
