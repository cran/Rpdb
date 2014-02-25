\name{range.coords}
\alias{range.coords}
\alias{range.atoms}
\alias{range.pdb}
\title{Range of Atomic Coordinates}
\description{
Determines the range of atomic coordinates.
}
\usage{
\method{range}{coords}(x, na.rm = FALSE, finite = FALSE, ...)

\method{range}{atoms}(x, na.rm = FALSE, finite = FALSE, ...)

\method{range}{pdb}(x, na.rm = FALSE, finite = FALSE, ...)
}
\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be omitted.}
  \item{finite}{logical, indicating if all non-finite elements should be omitted.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
Return a \code{\link{data.frame}} whose columns contain the range of the first, second and third coordinates of \code{x}.
}

\seealso{
\code{range}, \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
range(x)
range(range(x))

}
\keyword{ arith }
