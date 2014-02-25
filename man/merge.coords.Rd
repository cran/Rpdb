\name{merge.coords}
\alias{merge.coords}
\alias{merge.atoms}
\title{
Merge Two Objects of Class \sQuote{coords}
}
\description{
Merge by row two data frames of class \sQuote{coords}.
}
\usage{
\method{merge}{coords}(x, y, \dots)

\method{merge}{atoms}(x, y, reindex = TRUE, \dots)
}
\arguments{
  \item{x, y}{objects of class 'coords' to be merged.}
  \item{reindex}{a single element logical vector indicating if residue and element IDs have to be reindexed after merging.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{x} and \code{y} objects are merged together by row.
The \code{basis} attribute of these two objects must be the same (see \code{\link{basis}}).
For objects of class \sQuote{atoms} the residue and element IDs of \code{y}
are shifted to avoid any confusion with those of \code{x}.
}
\value{
Return a data.frame of the same class as \code{x} and \code{y} containing by row first \code{x} and then \code{y}.
If \code{x} and \code{y} have different \code{basis} attributes an error is returned.
}

\seealso{
\code{\link{coords}}, \code{\link{atoms}}, \code{\link{basis}}, \code{merge}, \code{merge.data.frame}, \code{\link{merge.pdb}}
}
\examples{
c1 <- coords( 1:3 ,  4:6 ,  7:9 , basis = "xyz")
c2 <- coords(10:12, 13:15, 16:18, basis = "xyz")
merge(c1,c2)

\dontrun{
## Merging objects with different basis sets return an error.
c2 <- coords(9:11, 12:14, 15:17, basis = "abc")
merge(c1,c2)
}

}
\keyword{ manip }
