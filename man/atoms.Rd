\encoding{UTF-8}
\name{atoms}
\alias{atoms}
\alias{atoms.default}
\title{
Create \sQuote{atoms} Object
}
\description{
Creates an object of class \sQuote{atoms} containing the data related to ATOM and HETATM records of a PDB file.
}
\usage{
atoms(...)

\method{atoms}{default}(recname, eleid, elename, alt, resname, chainid,
resid, insert, x1, x2, x3, occ, temp, segid, basis = "xyz", ...)
}
\arguments{
  \item{\dots}{arguments passed to methods.}
  \item{recname}{a character vector containing the record name for each element.}
  \item{eleid}{a integer vector containing the element ID for each element.}
  \item{elename}{a character vector containing the element name for each element.}
  \item{alt}{a character vector containing the alternate location indicator for each element.}
  \item{resname}{a character vector containing the residue name for each element.}
  \item{chainid}{a character vector containing the chain ID for each element.}
  \item{resid}{a integer vector containing the residue ID for each element.}
  \item{insert}{a character vector containing the codes for insertion of residue of each element.}
  \item{x1, x2, x3}{a numeric vector containing the first, second and third coordinate for each element.}
  \item{occ}{a numeric vector containing the occupancie for each element.}
  \item{temp}{a numeric vector containing the temperature factor for each element.}
  \item{segid}{a character vector containing the segment ID for each element.}
  \item{basis}{a single element character vector indicating the type of basis vector used to express the atomic coordinates.}
}
\details{
\code{atoms} is a generic function to create objects of class \sQuote{atoms}.
The purpose of this class is to store ATOM and HETATM records from PDB files.
The default method creates a \code{atoms} object from its different components, i.e.:
\code{recname}, \code{eleid}, \code{elename}, \code{alt}, \code{resname},
\code{chainid}, \code{resid}, \code{insert}, \code{x1}, \code{x2}, \code{x3},
\code{occ}, \code{temp}, \code{segid} and \code{basis}.
All the arguments have to be specified except \code{basis}
which by default is set to "xyz" (Cartesian coordinates).
}
\value{
Return a data.frame of class \sQuote{atoms} with the following components:
  \item{recname}{a character vector containing the record name for each element.}
  \item{eleid}{a integer vector containing the element ID for each element.}
  \item{elename}{a character vector containing the element name for each element.}
  \item{alt}{a character vector containing the alternate location indicator for each element.}
  \item{resname}{a character vector containing the residue name for each element.}
  \item{chainid}{a character vector containing the chain ID for each element.}
  \item{resid}{a integer vector containing the residue ID for each element.}
  \item{insert}{a character vector containing the codes for insertion of residue for each element.}
  \item{x1, x2, x3}{a numeric vector containing the first, second and third coordinate for each element.}
  \item{occ}{a numeric vector containing the occupencie for each element.}
  \item{temp}{a numeric vector containing the temperature factor for each element.}
  \item{segid}{a character vector containing the segment ID for each element.}
  \item{basis}{a single element character vector indicating the type of basis vector used to express the atomic coordinates.}
}

\seealso{
\code{\link{basis}}, \code{\link{coords}}, \code{\link{pdb}}
}
\examples{
x <- atoms(recname = c("ATOM","ATOM"), eleid = 1:2, elename = c("H","H"), alt = "",
            resname = c("H2","H2"), chainid = "", resid = c(1,1), insert = "",
            x1 = c(0,0), x2 = c(0,0), x3 = c(0,1), occ = c(0.0,0.0), temp = c(1.0,1.0),
            segid = c("H2","H2"))
print(x)
class(x)

}
\keyword{ classes }
