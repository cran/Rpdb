#' Number of Atoms in an Object Containing Atomic Coordinates
#' 
#' Evaluates the number of atoms in an object containing atomic coordinates.
#' 
#' \code{natom} is a generic function to evalute the number of atom in an object
#' containing atomic coordinates. The atomic coordinates of the object are first
#' filtered to keep ATOM and/or HETATM records as indicated by the 'ATOM' and 
#' 'HETATM' arguments. Then, if \code{factor} is specify, the object is splitted
#' to evalute the number of atoms in each group defined by \code{factor}. If 
#' \code{factor} is not specify then the total number of atoms in the object is 
#' return.
#' 
#' @return Return an integer or a vector of integer of lenght equal to 
#'   \code{nlevels(factor)} (if \code{factor} is specify) indication the number
#'   of atoms in the object or in the groups defined by \code{factor}.
#'   
#' @param x an R object containing atomic coordinates.
#' @param factor a factor used to split the object and evaluate the number of atoms
#'   in each group.
#' @param ATOM a single element logical vector indicating if ATOM records have
#'   to be considered or not.
#' @param HETATM a single element logical vector indicating if HETATM records
#'   have to be considered or not.
#' @param \dots further arguments passed to or from other methods.
#'   
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}},
#'   \code{\link[base]{factor}}, \code{\link[base]{split}}
#'   
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#'
#' natom(x)
#' natom(x, x$atoms$resid)
#' natom(x, x$atoms$resname)
#' natom(x, HETATM = FALSE)
#'   
#' @keywords manip
#' 
#' @name natom
#' @export
natom <- function(x, ...)
  UseMethod("natom")

#' @rdname natom
#' @export
natom.coords <- function(x, factor = NULL, ...)
{
  if(length(factor) == 0)
    nrow(x)
  else
    unlist(lapply(split(x, factor), nrow))
}

#' @rdname natom
#' @export
natom.atoms <- function(x, factor = NULL, ATOM = TRUE, HETATM = TRUE, ...)
{
  M <- ATOM & x$recname == "ATOM" | HETATM & x$recname == "HETATM"
  x <- x[M,]
  if(length(factor) != 0)
    factor <- factor[M]
  natom.coords(x, factor)
}

#' @rdname natom
#' @export
natom.pdb <- function(x, factor = NULL, ATOM = TRUE, HETATM = TRUE, ...)
  natom.atoms(x$atoms, factor, ATOM, HETATM, ...)
