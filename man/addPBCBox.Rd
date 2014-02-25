\name{addPBCBox}
\alias{addPBCBox}
\title{Add a PBC Box to the \sQuote{rgl} Scene}
\description{Add a PBC box to the current \sQuote{rgl} scene.}
\usage{
addPBCBox(x, lwd = 2)
}
\arguments{
  \item{x}{an object of class \sQuote{cryst1} containing unit cell parameters.}
  \item{lwd}{a numeric value indicating the line width used to draw the PBC box.}
}
\details{Add a box representing the Periodic Boundary Conditions of a molecular system.}
\value{
Return (using invisible) a two-column data.frame containing the IDs
and type indicators of the objects added to the scene.
}
\seealso{
\code{\link{addABC}}, \code{\link{addXYZ}}, \code{\link{visualize}}, \code{\link{rgl.cur}}, \code{\link{par3d}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
visualize(x, type = "l", abc = FALSE, pbc.box = FALSE)
addPBCBox(x$cryst1)
}
\keyword{ dynamic }
