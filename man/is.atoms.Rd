\name{is.atoms}
\alias{is.atoms}
\title{
Is an Object of Class \sQuote{atoms}?
}
\description{
Checks whether its argument is an object of class \sQuote{atoms}.
}
\usage{
is.atoms(x)
}
\arguments{
  \item{x}{
an R object to be tested
}
}
\details{
Tests if x is an object of class \sQuote{atoms}, i.e. if x has a \dQuote{class} attribute equal to \code{atoms}.
}
\value{
TRUE if x is an object of class \sQuote{atoms} and FALSE otherwise
}

\seealso{
\code{\link{atoms}}, \code{\link{pdb}}
}
\examples{
x <- atoms(recname = c("ATOM","ATOM"), eleid = 1:2, elename = c("H","H"), alt = "",
            resname = c("H2","H2"), chainid = "", resid = c(1,1), insert = "",
            x1 = c(0,0), x2 = c(0,0), x3 = c(0,1), occ = c(0.0,0.0), temp = c(1.0,1.0),
            segid = c("H2","H2"))
is.atoms(x)
is.atoms(1)

# Notice that the 'atoms' class extends the 'coords' class.
is.coords(x)

}
\keyword{ classes }
