#' From Cartesian to Fractional Coordinates and Vice Versa
#' 
#' Converts Cartesian coordinates into fractional coordinates and vice versa.
#' 
#' For \code{\link{atoms}} and \code{\link{pdb}} objects, the atomic coordinates
#' are first extracted from \code{x} using the \code{\link{coords}} function. 
#' Then, using the periodic boundary conditions stored into \code{crystal}, the 
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
#' @param crystal an object of class \code{\link{crystal}}.
#' @param cryst1 will be deprecated; use \code{crystal} instead.
#'   
#' @seealso \code{\link{basis}}, \code{\link{coords}}, \code{\link{atoms}},
#' \code{\link{pdb}}, \code{\link{crystal}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' basis(x)
#' x <- xyz2abc(x)
#' basis(x)
#' x <- abc2xyz(x)
#' basis(x)
#' 
#' \donttest{
#'   
#'   # This example returns an error because the coordinates stored
#'   # into the PDB file are already Cartesian coordinates.
#'   x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#'   try(x <- abc2xyz(x))
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
xyz2abc.coords <- function(x, crystal, ..., cryst1 = NULL)
{
	crystal = checkArgCrystal(crystal, cryst1);
	if(missing(crystal)) stop("Please specify a 'crystal' object");
	if(! is.coords(x)) stop("'x' must be an object of class 'coords'");
	if(! is.crystal(crystal))
		stop("'crystal' must be an object of class 'crystal'");
	
	if(basis(x) != "xyz") stop("Coordinates are not Cartesian coordinates")
	
	cell = cell.coords(crystal);
	to.return = solve(cell) %*% t(x);
	to.return = coords.default(
		to.return["a",], to.return["b",], to.return["c",], "abc");
	
	return(to.return);
}

#' @rdname xyz2abc
#' @export
xyz2abc.atoms <- function(x, crystal, ...)
{
	if(! is.atoms(x)) stop("'x' must be an object of class 'atoms'");
	
	value = xyz2abc.coords(coords.atoms(x), crystal);
	coords(x) = value;
	return(x);
}

#' @rdname xyz2abc
#' @export
xyz2abc.pdb <- function(x, crystal = x$crystal, ..., cryst1 = NULL)
{
	crystal = checkArgCrystal(crystal, cryst1);
	if(! is.pdb(x)) stop("'x' must be an object of class 'pdb'");
	value = xyz2abc.coords(coords.pdb(x), crystal);
	coords(x) = value;
	return(x);
}

#' @rdname xyz2abc
#' @export
xyz2abc.distances <- function(x, crystal, ...){
  if(!is.distances(x)) stop("'x' must be an object of class 'distances'")
  if(basis(x) != "xyz") stop("Coordinates are not Cartesian coordinates")
  
  dims <- dim(x$dx1)
  x <- coords(c(x$dx1), c(x$dx2), c(x$dx3), basis = basis(x))
  x <- xyz2abc.coords(x, crystal = crystal)
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
abc2xyz.coords <- function(x, crystal, ..., cryst1 = NULL)
{
	crystal = checkArgCrystal(crystal, cryst1);
	if(missing(crystal)) stop("Please specify a 'crystal' object");
	if(! is.coords(x)) stop("'x' must be an object of class 'coords'");
	if(! is.crystal(crystal)) stop("'crystal' must be an object of class 'crystal'");
	
	if(basis.default(x) != "abc")
		stop("Coordinates are not fractional coordinates");
	
	cell <- cell.coords(crystal);
	to.return <- cell %*% t(x);
	to.return <- coords.default(
		to.return["x",], to.return["y",], to.return["z",], "xyz");
  
	return(to.return)
}

#' @rdname xyz2abc
#' @export
abc2xyz.atoms <- function(x, crystal, ...)
{
	if(! is.atoms(x)) stop("'x' must be an object of class 'atoms'")
	
	value = abc2xyz.coords(coords.atoms(x), crystal);
	coords(x) = value;
	return(x)
}

#' @rdname xyz2abc
#' @export
abc2xyz.pdb <- function(x, crystal = x$crystal, ..., cryst1 = NULL)
{
	crystal = checkArgCrystal(crystal, cryst1);
	if(! is.pdb(x)) stop("'x' must be an object of class 'pdb'");
	
	value = abc2xyz.coords(coords.pdb(x), crystal);
	coords(x) = value;
	return(x)
}

#' @rdname xyz2abc
#' @export
abc2xyz.distances <- function(x, crystal, ...){
	if(! is.distances(x)) stop("'x' must be an object of class 'distances'")
	if(basis(x) != "abc") stop("Coordinates are not fractional coordinates")
	
	dims = dim(x$dx1);
	x   = coords(c(x$dx1), c(x$dx2), c(x$dx3), basis = basis(x));
	x   = abc2xyz.coords(x, crystal = crystal);
	dx1 = array(x$x1, dim = dims);
	dx2 = array(x$x2, dim = dims);
	dx3 = array(x$x3, dim = dims);
	x   = distances.default(dx1, dx2, dx3, basis = "xyz");
	return(x)
}
