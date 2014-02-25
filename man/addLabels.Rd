\name{addLabels}
\alias{addLabels}
\alias{info3d}
\alias{info3d.atoms}
\alias{info3d.pdb}
\alias{addResLab}
\alias{addResLab.atoms}
\alias{addResLab.pdb}
\alias{addEleLab}
\alias{addEleLab.atoms}
\alias{addEleLab.pdb}
\title{Add Labels to the \sQuote{rgl} Scene}
\description{Add Labels to the current \sQuote{rgl} scene.}
\usage{
addResLab(x, ...)
\method{addResLab}{atoms}(x, at.centre = TRUE, col = "black", ...)
\method{addResLab}{pdb}(x, at.centre = TRUE, col = "black", ...)

addResLab(x, ...)
\method{addEleLab}{atoms}(x, eleid = FALSE, col = "black", ...)
\method{addEleLab}{pdb}(x, eleid = FALSE, col = "black", ...)

info3d(...)
\method{info3d}{atoms}(x, id = rgl.ids(), col = "black", verbose = TRUE, adj = 0, ...)
\method{info3d}{pdb}(x, id = rgl.ids(), col = "black", verbose = TRUE, adj = 0, ...)
}
\arguments{
  \item{x}{an R object containing atomic coordinates.}
  \item{at.centre}{a single element logical vector indicating if residue labels have to be added only at the position of the residue's centre-of-mass instead of at each atomic position.}
  \item{eleid}{a single element logical vector indicating if the element ids have to be concatenated with the element names to prepare the labels.}
  \item{col}{the colors used to display the labels.}
  \item{adj}{one value specifying the horizontal adjustment, or two, specifying horizontal and vertical adjustment respectively. See \code{\link{rgl.texts}}}
  \item{id}{vector of ID numbers of \sQuote{rgl} items, as returned by \code{rgl.ids}. The vertexes of these items are used to display the labels.}
  \item{verbose}{a logical value specifying if information have to be printed to the terminal.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{addResLab} add residue labels to the scene. If \code{at.centre==TRUE} only one label per residue is added at the centre of the residue. Otherwise, residue labels are added at each atomic positions.

\code{addEleLab} add element labels to the scene at each atomic positions.

\code{info3d} activate an interactive mode to add labels by selecting atoms by \bold{right-clicing} on the current \sQuote{rgl} scene. To escape the active mode press the ESC key.
}
\value{
\code{addResLab} and \code{addEleLab} return (using invisible) a two-column data.frame containing the IDs
and type indicators of the objects added to the scene.
}

\seealso{
\code{\link{pdb}}, \code{\link{visualize}}, \code{\link{measure}}
}
\examples{
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
visualize(x, type = "l", mode = NULL)
addResLab(x)
x <- read.pdb(system.file("examples/Pentacene.pdb",package="Rpdb"))
visualize(x, type = "l", mode = NULL)
addEleLab(x)
}
\keyword{ dynamic }
