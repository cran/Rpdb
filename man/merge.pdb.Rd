\name{merge.pdb}
\alias{merge.pdb}
\title{
Merge Two Obejct of Class \sQuote{pdb}
}
\description{
Merge two objects of class \sQuote{pdb}.
}
\usage{
\method{merge}{pdb}(x, y, reindex = TRUE, ...)
}
\arguments{
  \item{x, y}{objects of class \sQuote{pdb} to be merged.}
  \item{reindex}{a single element logical vector indicating if residue and element IDs should be reindexed after merging.}
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
The \code{atoms} and \code{conect} components of the two \code{pdb} objects are merged together by row.
The \code{basis} attribute of \code{x} and \code{y} must be the same (see \code{\link{basis}}).
The residue and element IDs of \code{y} are shifted to avoid any confusion with those of \code{x}.
If the \code{cryst1} components of \code{x} and \code{y} differ the \code{cryst1} component of \code{x}
is used to build the returned object.
}
\value{
Return an object of class \sQuote{pdb} merging \code{x} and \code{y}.
If \code{x} and \code{y} have different \code{basis} attributes an error is returned.
}

\seealso{
\code{\link{pdb}}, \code{\link{basis}}, \code{merge}, \code{\link{merge.coords}}
}
\examples{
## Merge two 'pdb' objects
x1 <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
x2 <- x1
x2$atoms$x1 <- x2$atoms$x1 + 10
y <- merge(x1, x2)

## Two merge more than two 'pdb' objects together
## first create a list of 'pdb' objects and use the reduce function to merged them all.
x3 <- x1
x3$atoms$x1 <- x3$atoms$x1 - 10
L <- list(x1, x2, x3)
y <- Reduce(merge, L)

}
\keyword{ manip }
