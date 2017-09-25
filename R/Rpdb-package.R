#' Read, Write, Visualize and Manipulate PDB Files
#' 
#' Provides tools to read, write, visualize PDB files, and perform structural manipulations.
#' 
#' This package has been developed for computational chemists wishing
#' to manipulate molecular structures stored in PDB files. PDB files can easily be
#' read, written, visualized and some basic structural manipulations can be achieved with the
#' present package.
#' Convertion of Cartesian coordinates into fractional coordinates.
#' Spliting a molecular structure into fragments.
#' Computation of centers-of-geometry and centers-of-mass.
#' Wrapping molecular structure using periodical boundary conditions.
#' Tranlation, rotation and reflection of atomic coordinates.
#' Calculate atomic bond lengths, angles and dihedrals. 
#'
#' @author Julien Idé \email{julien.ide.fr@gmail.com}\cr
#'
#' @references
#' More information on the PDB format can be found here:\cr
#' http://www.wwpdb.org/documentation/format33/v3.3.html
#'
#' @keywords package
#' 
#' @examples
#' ## Read a PDB file included in the package
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#' 
#' ## Visualize the PDB file
#' visualize(x, mode = NULL)
#' 
#' ## From Cartesian to fractional coordinates and vice versa
#' x <- xyz2abc(x)
#' basis(x)
#' natom(x,x$atoms$resid)
#' range(x)
#' centres(x)
#' x <- abc2xyz(x)
#' basis(x)
#' natom(x,x$atoms$resid)
#' range(x)
#' centres(x)
#' 
#' ## Split and unsplit
#' F <- x$atoms$resid
#' x <-   split(x, F)
#' x <- unsplit(x, F)
#' 
#' ## Subset and merge
#' x.PCB.only <- subset(x, resname == "PCB")
#' x.DCB.only <- subset(x, resname == "DCB")
#' x <- merge(x.PCB.only, x.DCB.only)
#' 
#' ## Duplicate and wrap
#' x <- replicate(x, a.ind = -1:1, b.ind = -1:1, c.ind = -1:1)
#' x <- wrap(x)
#' 
#' ## Write the 'pdb' object 'x' in a temporary file.
#' write.pdb(x, file = tempfile())
#' 
#' @name Rpdb-package
NULL
