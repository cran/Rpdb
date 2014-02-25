\name{addXYZ}
\alias{addXYZ}
\title{Add Cartesian Axes to the \sQuote{rgl} Scene}
\description{Add Cartesian axes to the current \sQuote{rgl} scene.}
\usage{
addXYZ(lwd = 2, labels = TRUE, cex = 2)
}
\arguments{
  \item{lwd}{a numeric value indicating the line width used to draw the axes.}
  \item{labels}{a logical value indicating whether the labels of the axes have to bo drawn.}
  \item{cex}{a numeric value indicating the magnification used to draw the labels of the axes.}
}
\details{Add the Cartesian axes x, y and z to the current rgl device.}

\value{
Return (using invisible) a two-column data.frame containing the IDs
and type indicators of the objects added to the scene.
}
\seealso{
\code{\link{addABC}}, \code{\link{addPBCBox}}, \code{\link{visualize}}, \code{\link{rgl.cur}}, \code{\link{par3d}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
visualize(x, type = "l")
addXYZ()
}
\keyword{ dynamic }
