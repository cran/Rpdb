#' From Cartesian to Fractional Coordinates and Vis Versa
#' 
#' Converts Cartesian coordinates into fractional coordinates and vice versa.
#' 
#' For \code{\link{atoms}} and \code{\link{pdb}} objects, the atomic coordinates
#' are first extracted from \code{x} using the \code{\link{coords}} function. 
#' Then, using the periodic boundary conditions stored into \code{cryst1}, the 
#' coordinates are converted from Cartesian to fractional (for the 
#' \code{xyz2abc} functions) or from fractional to Cartesian (for the 
#' \code{abc2xyz} functions) coordinates. Finally, for \code{\link{atoms}} and 
#' \code{\link{pdb}} objects, the new atomic coordinates are reassigned to the 
#' original \code{x} object using the \code{\link{coords<-}} function and 
#' \code{x} is returned.
#' 
#' @return Return an object of the same class as \code{x}, with atomic 
#'   coordinates expressed in a different basis set.
#'   
#' @param \dots arguments passed to methods.
#' @param x an R object containing atomic coordinates.
#' @param cryst1 an object of class \code{\link{cryst1}}.
#'   
#' @seealso \code{\link{basis}}, \code{\link{coords}}, \code{\link{atoms}},
#' \code{\link{pdb}}, \code{\link{cryst1}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#' basis(x)
#' x <- xyz2abc(x)
#' basis(x)
#' x <- abc2xyz(x)
#' basis(x)
#' 
#' \dontrun{
#'   
#'   # This example return an error because the coordinates stored
#'   # into the PDB file are already Cartesian coordinates.
#'   x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#'   x <- abc2xyz(x)
#' }
#' 
#' @keywords manip
#' 
#' @name xyz2abc
#' @export
xyz2abc <- function(...)
  UseMethod("xyz2abc")

#' @rdname xyz2abc
#' @export
xyz2abc.coords <- function(x, cryst1, ...)
{
  if(missing(cryst1)) stop("Please specify a 'cryst1' object")
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  if(basis(x) != "xyz") stop("Coordinates are not Cartesian coordinates")
  
  cell <- cell.coords(cryst1)
  to.return <- solve(cell)%*%t(x)
  to.return <- coords.default(to.return["a",],to.return["b",],to.return["c",],"abc")

  return(to.return)
  
}

#' @rdname xyz2abc
#' @export
xyz2abc.atoms <- function(x, cryst1, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  value <- xyz2abc.coords(coords.atoms(x), cryst1)
  coords(x) <- value
  
  return(x)
}

#' @rdname xyz2abc
#' @export
xyz2abc.pdb <- function(x, cryst1 = x$cryst1, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  value <- xyz2abc.coords(coords.pdb(x), cryst1)
  coords(x) <- value
  
  return(x)
}

#' @rdname xyz2abc
#' @export
xyz2abc.distances <- function(x, cryst1, ...){
  if(!is.distances(x)) stop("'x' must be an object of class 'distances'")
  if(basis(x) != "xyz") stop("Coordinates are not Cartesian coordinates")
  
  dims <- dim(x$dx1)
  x <- coords(c(x$dx1), c(x$dx2), c(x$dx3), basis = basis(x))
  x <- xyz2abc.coords(x, cryst1 = cryst1)
  dx1 <- array(x$x1, dim = dims)
  dx2 <- array(x$x2, dim = dims)
  dx3 <- array(x$x3, dim = dims)
  x <- distances.default(dx1, dx2, dx3, basis = "abc")
  return(x)
}

#' @rdname xyz2abc
#' @export
abc2xyz <- function(...)
  UseMethod("abc2xyz")

#' @rdname xyz2abc
#' @export
abc2xyz.coords <- function(x, cryst1, ...)
{
  if(missing(cryst1)) stop("Please specify a 'cryst1' object")
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  if(basis.default(x) != "abc") stop("Coordinates are not fractional coordinates")
  
  cell <- cell.coords(cryst1)
  to.return <- cell%*%t(x)
  to.return <- coords.default(to.return["x",],to.return["y",],to.return["z",],"xyz")
  
  return(to.return)
  
}

#' @rdname xyz2abc
#' @export
abc2xyz.atoms <- function(x, cryst1, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  value <- abc2xyz.coords(coords.atoms(x), cryst1)
  coords(x) <- value
  
  return(x)
}

#' @rdname xyz2abc
#' @export
abc2xyz.pdb <- function(x, cryst1 = x$cryst1, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  value <- abc2xyz.coords(coords.pdb(x), cryst1)
  coords(x) <- value
  return(x)
}

#' @rdname xyz2abc
#' @export
abc2xyz.distances <- function(x, cryst1, ...){
  if(!is.distances(x)) stop("'x' must be an object of class 'distances'")
  if(basis(x) != "abc") stop("Coordinates are not fractional coordinates")
  
  dims <- dim(x$dx1)
  x <- coords(c(x$dx1), c(x$dx2), c(x$dx3), basis = basis(x))
  x <- abc2xyz.coords(x, cryst1 = cryst1)
  dx1 <- array(x$x1, dim = dims)
  dx2 <- array(x$x2, dim = dims)
  dx3 <- array(x$x3, dim = dims)
  x <- distances.default(dx1, dx2, dx3, basis = "xyz")
  return(x)
}
