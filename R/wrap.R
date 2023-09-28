#' Wrap Atomic Coordinates
#' 
#' Wraps atomic coordinates using periodic boundary conditions.
#' 
#' The \code{wrap} function translates all atoms out of the unit cell back into
#' the unit cell using periodic boundary conditions. To do so, the \code{wrap}
#' function first converts Cartesian into fractional coordinates. Then atoms
#' with fractional coordinates greater than 1 or lower than 0 are respectively
#' translated by -1 or +1. Finally, if the original atomic coordinates were
#' Cartesian coordinates their are reconverted into Cartesian coordinates.
#' 
#' @return Return an object of class \sQuote{pdb} with wrapped atomic
#'   coordinates.
#'   
#' @param x an R object containing atomic coordinates to be wrapped.
#' @param cryst1 an object of class \sQuote{cryst1} containing periodic boundary
#'   conditions used for wrapping.
#' @param factor a factor used to wrap the atoms by groups.
#' @param \dots further arguments passed to or from other methods.
#'   
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}},
#' \code{\link{cryst1}}, \code{\link{centres.pdb}}, \code{\link{xyz2abc}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' 
#' # Translation of the atoms along x-axis
#' x$atoms$x1 <- x$atoms$x1 + 10
#' 
#' # Wrapping the structure
#' y <- wrap(x)
#' 
#' @keywords manip
#' 
#' @name wrap
#' @export
wrap <- function(x, ...)
  UseMethod("wrap")

#' @rdname wrap
#' @export
wrap.coords <- function(x, cryst1 = NULL, factor = NULL, ...)
{
  if(is.null(cryst1)) stop("Please specify a 'cryst1' object")
  if(is.null(factor)) factor <- 1:natom(x)
  
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")

  b <- basis(x)
  if(b == "xyz") x <- xyz2abc(x, cryst1)
  
  centers <- centres.coords(x, factor = factor, unsplit = TRUE)
  x[centers > 1] <- x[centers > 1] - 1
  x[centers < 0] <- x[centers < 0] + 1
  
  if(b == "xyz") x <- abc2xyz.coords(x, cryst1)
  
  return(x)
}

#' @rdname wrap
#' @export
wrap.atoms <- function(x, cryst1= NULL, factor = NULL, ...)
{
  if(is.null(cryst1)) stop("Please specify a 'cryst1' object")
  if(is.null(factor)) factor <- x$resid
  
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  coords(x) <- wrap.coords(coords(x), cryst1, factor)
  
  return(x)
}

#' @rdname wrap
#' @export
wrap.pdb <- function(x, cryst1 = x$cryst1, factor = NULL, ...)
{
  if(is.null(factor)) factor <- x$atoms$resid
  
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  coords(x) <- wrap.coords(coords(x), cryst1, factor)
  
  return(x)
}
