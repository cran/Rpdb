% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.pdb.R
\name{read.pdb}
\alias{read.pdb}
\title{PDB File Reader}
\usage{
read.pdb(
  file,
  ATOM = TRUE,
  HETATM = TRUE,
  CRYST1 = TRUE,
  CONECT = TRUE,
  TITLE = TRUE,
  REMARK = TRUE,
  MODEL = 1
)
}
\arguments{
\item{file}{a single element character vector containing the name of the PDB file to be read.}

\item{ATOM}{a single element logical vector indicating whether ATOM records have to be read.}

\item{HETATM}{a single element logical vector indicating whether HETATM records have to be read.}

\item{CRYST1}{a single element logical vector indicating whether CRYST1 records have to be read.}

\item{CONECT}{a single element logical vector indicating whether CONECT records have to be read.}

\item{TITLE}{a single element logical vector indicating whether TITLE records have to be read.}

\item{REMARK}{a single element logical vector indicating whether REMARK records have to be read.}

\item{MODEL}{an integer vector containing the serial number of the MODEL sections to be read. Can also be equal to NULL to read all the MODEL sections or to NA to ignore MODEL records (see details).}
}
\value{
When a single MODEL section is read, this function returns an object of class  \sQuote{pdb} (a list with a \code{class} attribute equal to \code{pdb}) with the following components:
\item{title}{a character vector containing the TITLE records found in the PDB file.}
\item{remark}{a character vector containing the REMARK records found in the PDB file.}
\item{cryst1}{a list of class \sQuote{cryst1} containing the first CRYST1 record found in the PDB file. All others are ignored.}
\item{atoms}{a data.frame of class \sQuote{atoms} containing the ATOM and HETATM records found in the PDB file.}
\item{conect}{a data.frame of class \sQuote{conect} containing the CONECT records found in the PDB file.}
When multiple MODEL sections are read, a list of object of class \sQuote{pdb} is returned.
}
\description{
Reads a Protein Data Bank (PDB) coordinate file.
}
\details{
The \code{read.pdb} function reads the TITLE, REMARK, ATOM, HETATM, CRYST1 and CONECT records from a PDB file. Three different reading modes can be used depending on the value of \code{MODEL}: 
\itemize{
  \item When \code{MODEL} is a vector of integers, MODEL sections whose serial numbers match these integers are read.
  \item When \code{MODEL == NULL}, all MODEL sections are read.
  \item When \code{MODEL == NA}, MODEL records are ignored to read all ATOM and/or HETATM records together to return a single object.
}
}
\examples{
# Read a PDB file included with the package
x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))

# Visualize the PDB file
visualize(x, mode = NULL)

# Write the 'pdb' object 'x' in file "Rpdb.pdb" into the current directory
write.pdb(x, file = "Rpdb.pdb")

# Cleanup
unlink("Rpdb.pdb")

}
\references{
PDB format has been taken from:
http://www.wwpdb.org/documentation/format33/v3.3.html
}
\seealso{
\code{\link{write.pdb}}, \code{\link{pdb}}, \code{\link{cryst1}}, \code{\link{atoms}}, \code{\link{conect}}
}
\keyword{IO}
