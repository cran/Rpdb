#' Rotation of Atomic Coordinates
#' 
#' Rotation of atomic coordinates around a given vector.
#' 
#' \code{R} is generic functions. Method for objects of class \sQuote{coords}
#' first convert the coordinates into Cartesian coordinates using \code{cryst1}
#' if needed. Once rotated, the coordinates are reconverted back to the orginal
#' basis set using again \code{cryst1}. Method for objects of class \sQuote{pdb}
#' first extract coordinates from the object using the function \code{coords},
#' perform the rotation, and update the coordinates of the \sQuote{pdb} object
#' using the function \code{coords<-}.
#' 
#' @return An object of the same class as \code{x} with rotated coordinates.
#' 
#' @param obj an R object containing atomic coordinates.
#' @param angle the angle of the rotation in degrees.
#' @param x the x-component of the rotation vector.
#' @param y the y-component of the rotation vector.
#' @param z the z-component of the rotation vector.
#' @param mask a logical vector indicating the set of coordinates to which the rotation has to be applyed.
#' @param cryst1 an object of class \sQuote{cryst1} use to convert fractional into Cartesian coordinates when need.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso 
#' Helper functions for rotation around a given Cartesian vector:\cr
#' \code{\link{Rx}}, \code{\link{Ry}}, \code{\link{Rz}}\cr
#' Passing from Cartesian to fractional coordinates (or Vis Versa):\cr
#' \code{\link{xyz2abc}}, \code{\link{abc2xyz}}
#' 
#' @examples 
#' # First lets read a pdb file
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#' cell <- cell.coords(x)
#' visualize(x, mode = NULL)
#' # Rotation of the structure around the c-axis
#' visualize(R(x, 90, x=cell["x","c"], y=cell["y","c"], z=cell["z","c"]),
#'           mode = NULL)
#' # Rotation of the residue 1 around the c-axis
#' visualize(R(x, 90, x=cell["x","c"], y=cell["y","c"], z=cell["z","c"], mask=x$atoms$resid==1),
#'           mode = NULL)
#'           
#' @keywords manip
#' 
#' @name rotation
#' @export
R <- function(...)
  UseMethod("R")

#' @rdname rotation
#' @export
R.coords <- function(obj, angle = 0, x = 0, y = 0, z = 1, mask = TRUE, cryst1 = NULL, ...){
  if(!is.coords(obj)) stop("'object' must be an obj of class 'coords'")

  if(length(mask) != natom(obj)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(obj))
  }
  
  basis.ori <- basis(obj)
  if(basis.ori != "xyz"){
    if(is.null(cryst1))
      stop("Please specify a 'cryst1' obj to convert your fractional into Cartesian coordinates")
    obj <- abc2xyz(obj, cryst1 = cryst1)
  }
  M <- rgl::rotationMatrix(angle=angle*pi/180, x = x, y = y, z = z)[1:3,1:3]
  obj[mask,] <- coords(as.matrix(obj[mask,])%*%M, basis = "xyz")
  if(basis.ori != "xyz")
    obj <- xyz2abc(obj, cryst1 = cryst1)
  return(obj)
}

#' @rdname rotation
#' @export
R.pdb <- function(obj, angle = 0, x = 0, y = 0, z = 1, mask = TRUE, cryst1 = obj$cryst1, ...){
  if(!is.pdb(obj)) stop("'object' must be an obj of class 'pdb'")
  
  coords(obj) <- R(coords(obj), angle = angle, x = x, y = y, z = z, mask = mask, cryst1 = cryst1, ...)
  
  return(obj)
}
