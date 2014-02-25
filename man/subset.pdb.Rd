\name{subset.pdb}
\alias{subset.pdb}
\alias{subset.atoms}
\title{
Subsetting \sQuote{atoms} and \sQuote{pdb} Objects
}
\description{
Return subsets of \sQuote{atoms} or \sQuote{pdb} objects which meet conditions.
}
\usage{
\method{subset}{pdb}(x, subset, drop = FALSE, reindex.all = TRUE, ...)

\method{subset}{atoms}(x, subset, drop = FALSE, reindex.all = TRUE, ...)
}
\arguments{
  \item{x}{object to be subsetted.}
  \item{subset}{logical expression indicating elements or rows to keep: missing values are taken as false.}
  \item{drop}{passed on to [ indexing operator.}
  \item{reindex.all}{a single element logical vector indicating if residues and elements IDs have to be reindexed after subsetting.}
  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
For a \sQuote{atoms} object the method is similar to the data.frame method (see \code{\link{subset}}) but allow to directly reindex the elements and residues IDs. For a \sQuote{pdb} object subsetting is apply on the \code{atoms} and \code{conect} components of the object in a consistent way. First the \code{atoms} component is subsetted and then the \code{conect} component is filtered to keep only the conectivity for the subset.  
}
\value{
Return a subsetted object of the same class as \code{x}.
}

\seealso{
\code{\link{subset}}, \code{\link{pdb}}, \code{\link{atoms}}, \code{\link{reindex}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
y <- subset(x, x$atoms$eleid \%in\% sample(x$atoms$eleid, 10))
is(y)
y <- subset(x$atoms, x$atoms$eleid \%in\% sample(x$atoms$eleid, 10))
is(y)
x <- coords(x)
y <- subset(x, x < 0)
is(y)

}
\keyword{ manip }
