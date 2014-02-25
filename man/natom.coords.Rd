\name{natom.coords}
\alias{natom.coords}
\alias{natom.atoms}
\alias{natom.pdb}
\title{
Number of Atoms in \sQuote{coords}, \sQuote{atoms} or \sQuote{pdb} Objects
}
\description{Evaluates the number of atoms in an object of class \sQuote{coords}, \sQuote{atoms} or \sQuote{pdb}.}
\usage{
\method{natom}{coords}(x, factor = NULL, \dots)

\method{natom}{atoms}(x, factor = NULL, ATOM = TRUE, HETATM = TRUE, \dots)

\method{natom}{pdb}(x, factor = NULL, ATOM = TRUE, HETATM = TRUE, \dots)
}

\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{factor}{a factor use to split the object and evalute the number of atom in each group.}
  \item{ATOM}{a single element logical vector indicating if ATOM records have to be considered or not.}
  \item{HETATM}{a single element logical vector indicating if HETATM records have to be considered or not.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
The atomic coordinates of the object are first filtered to keep ATOM and/or HETATM records as indicated by the 'ATOM' and 'HETATM' arguments. Then, if \code{factor} is specify, the object is splitted to evalute the number of atoms in each group defined by \code{factor}. If \code{factor} is not specify then the total number of atoms in the object is return.
}
\value{
Return an integer or a vector of integer of lenght equal to \code{nlevels(factor)} (if \code{factor} is specify) indication the number of atoms in the object or in the groups defined by \code{factor}.
}

\seealso{
\code{\link{natom}}, \code{\link{factor}}, \code{\link{split}}, \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}} 
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
natom(x)
natom(x, x$atoms$resid)
natom(x, x$atoms$resname)
natom(x, HETATM=FALSE)

}
\keyword{ array }
