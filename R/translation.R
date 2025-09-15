#' Translation of Atomic Coordinates
#' 
#' Translation of Cartesian or fractional coordinates.
#' 
#' \code{Txyz} and \code{Tabc} are generic functions. Method for objects of 
#' class \sQuote{coords} first convert the coordinates into Cartesian or 
#' fractional coordinates using \code{crystal} if needed to performed the 
#' translation. Once translated, the coordinates are reconverted back to the 
#' original basis set using again \code{crystal}. Method for objects of class 
#' \sQuote{pdb} first extract coordinates from the object using the function 
#' \code{coords}, perform the translation, and update the coordinates of the 
#' \sQuote{pdb} object using the function \code{coords<-}. The \code{thickness} 
#' argument can be use to translate selected atoms by a fraction of its 
#' thickness along the translation direction. This can be use when merging two 
#' fragments centered at the origin to build a dimer to avoid atomic overlap and
#' set the inter-fragment distance (see examples).
#' 
#' @return An object of the same class as \code{x} with translated coordinates.
#'   
#' @param obj an R object containing atomic coordinates.
#' @param x the x-component of the translation vector.
#' @param y the y-component of the translation vector.
#' @param z the z-component of the translation vector.
#' @param a the a-component of the translation vector.
#' @param b the b-component of the translation vector.
#' @param c the c-component of the translation vector.
#' @param mask a logical vector indicating the set of coordinates to which to
#'   apply the translation.
#' @param thickness a numeric value indicating the fraction of the thicknees of
#'   the selected atom to be added to the translation vector (Usually 0, 0.5 or
#'   1. See details).
#' @param cryst1 an object of class \sQuote{crystal} use to convert Cartesian
#'   into fraction coordinates (or Vis Versa) when need.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso Helper functions for translation along given Cartesian or lattice vector:\cr
#' \code{\link{Tx}}, \code{\link{Ty}}, \code{\link{Tz}}, \code{\link{Ta}}, \code{\link{Tb}}, \code{\link{Tc}}\cr
#' Passing from Cartesian to fractional coordinates (or Vis Versa):\cr
#' \code{\link{xyz2abc}}, \code{\link{abc2xyz}}
#' 
#' @examples 
#' # First lets read a pdb file
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' visualize(x, mode = NULL)
#' visualize(Txyz(x, y=10), mode = NULL)
#' visualize(Txyz(x, y=10, mask = x$atoms$resid==1), mode = NULL)
#' visualize(Tabc(x, b=1 ), mode = NULL)
#' visualize(Tabc(x, b=1 , mask = x$atoms$resid==1), mode = NULL)
#' 
#' # Lets build a C70/Pentacene dimer with an inter-molecular distance equal to 3.5
#' C70 <- read.pdb(system.file("examples/C70.pdb", package="Rpdb"))
#' Pen <- read.pdb(system.file("examples/Pentacene.pdb", package="Rpdb"))
#' x <- merge(C70, Pen)
#' visualize(x, mode = NULL)
#' viewXY()
#' visualize(Txyz(x, x=0, y=0, z=3.5, mask = x$atoms$resname == "C70", thickness=0.5), mode = NULL)
#' viewXY()
#' 
#' @keywords manip
#' 
#' @name translation
#' @export
Txyz <- function(...)
  UseMethod("Txyz")

#' @rdname translation
#' @export
Txyz.coords <- function(obj, x = 0, y = 0, z = 0, mask = TRUE,
		thickness = NULL, cryst1 = NULL, ...) {
  if(!is.coords(obj)) stop("'object' must be an obj of class 'coords'")
  
  if(length(mask) != natom(obj)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(obj))
  }
  
  v <- coords(x,y,z, basis = "xyz")
  T <- coords(0,0,0, basis = "xyz")
  if(basis(obj) != "xyz"){
    if(is.null(cryst1))
      stop("Please specify a 'crystal' obj to convert the fractional coordinates into Cartesian")
    v <- xyz2abc(v, crystal = cryst1)
    T <- xyz2abc(T, crystal = cryst1)
  }

  vn <- coords(0,0,0, basis = "xyz")
  if(sqrt(sum(v^2)) != 0) vn <- v/sqrt(sum(v^2))
  
  if(!is.null(thickness)) {
    if(length(thickness) != 1) stop("'thickness must be a single element numeric vector'")
    T <- as.matrix(obj[mask,])%*%t(vn)
    T <- diff(range(T))*vn*thickness
  }
  
  obj$x1[mask] <- obj$x1[mask] + v$x1 + T$x1
  obj$x2[mask] <- obj$x2[mask] + v$x2 + T$x2
  obj$x3[mask] <- obj$x3[mask] + v$x3 + T$x3

  return(obj)
}

#' @rdname translation
#' @export
Txyz.pdb <- function(obj, x = 0, y = 0, z = 0, mask = TRUE,
		thickness = NULL, cryst1 = obj$crystal, ...) {
  if(!is.pdb(obj)) stop("'object' must be an obj of class 'pdb'")
  
  coords(obj) <- Txyz(coords(obj), x = x, y = y, z = z,
		mask = mask, thickness = thickness, cryst1 = cryst1, ...);
  
  return(obj)
}

#' @rdname translation
#' @export
Tabc <- function(...)
  UseMethod("Tabc")

#' @rdname translation
#' @export
Tabc.coords <- function(obj, a = 0, b = 0, c = 0, mask = TRUE,
		thickness = NULL, cryst1 = NULL, ...){  
  if(!is.coords(obj)) stop("'object' must be an obj of class 'coords'")

  if(length(mask) != natom(obj)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(obj))
  }
  
  v <- coords(a,b,c, basis = "abc")
  T <- coords(0,0,0, basis = "abc")
	if(basis(obj) != "abc"){
		if(is.null(cryst1))
			stop("Please specify a 'crystal' obj to convert the Cartesian into fractional coordinates");
		v = abc2xyz(v, crystal = cryst1);
		T = abc2xyz(T, crystal = cryst1);
	}

  vn <- coords(0,0,0, basis = "abc")
  if(sqrt(sum(v^2)) != 0) vn <- v/sqrt(sum(v^2))

  if(!is.null(thickness)) {
    if(length(thickness) != 1) stop("'thickness must be a single element numeric vector'")
    T <- as.matrix(obj[mask,])%*%t(vn)
    T <- diff(range(T))*vn*thickness
  }
  
  obj$x1[mask] <- obj$x1[mask] + v$x1 + T$x1
  obj$x2[mask] <- obj$x2[mask] + v$x2 + T$x2
  obj$x3[mask] <- obj$x3[mask] + v$x3 + T$x3
  
  return(obj)
}

#' @rdname translation
#' @export
Tabc.pdb <- function(obj, a = 0, b = 0, c = 0, mask = TRUE,
		thickness = NULL, cryst1 = obj$crystal, ...) {
	if(! is.pdb(obj)) stop("'object' must be an obj of class 'pdb'");
	
	coords(obj) = Tabc(coords(obj), a = a, b = b, c = c,
		mask = mask, thickness = thickness, cryst1 = cryst1, ...);
	
	return(obj)
}
