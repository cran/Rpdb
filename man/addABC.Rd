\name{addABC}
\alias{addABC}
\title{Add Lattice Vectors to the \sQuote{rgl} Scene}
\description{Add lattice vectors to the current \sQuote{rgl} scene.}
\usage{
addABC(x, lwd = 2, labels = TRUE, cex = 2)
}
\arguments{
  \item{x}{an object of class \sQuote{cryst1} containing unit cell parameters.}
  \item{lwd}{a numeric value indicating the line width used to draw the lattice vectors.}
  \item{labels}{a logical value indicating whether the labels of the lattice vectors have to be drawn.}
  \item{cex}{a numeric value indicating the magnification used to draw the labels of the lattice vectors.}
}
\details{Add the lattice vectors a, b and c to the current rgl device.}
\value{
Return (using invisible) a two-column data.frame containing the IDs
and type indicators of the objects added to the scene.
}

\seealso{
\code{\link{addXYZ}}, \code{\link{addPBCBox}}, \code{\link{visualize}}, \code{\link{rgl.cur}}, \code{\link{par3d}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
visualize(x, type = "l", abc = FALSE, pbc.box = FALSE)
addABC(x$cryst1)
}
\keyword{ dynamic }
