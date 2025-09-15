#' Reflexion of Atomic Coordinates
#' 
#' Perform a reflexion (or mirror) operation on atomic coordinates with respect 
#' to a given reflexion plane.
#' 
#' \code{mirror} is a generic function.
#' Method for objects of class \sQuote{coords} first convert the coordinates
#' into Cartesian coordinates using \code{crystal} if needed.
#' Once reflected, the coordinates are reconverted back to the original basis set
#' using again \code{crystal}.
#'
#' Method for objects of class \sQuote{pdb} first extract coordinates from the object 
#' using the function \code{coords}, perform the reflection, and update the 
#' coordinates of the \sQuote{pdb} object using the function \code{coords<-}.
#' 
#' @return An object of the same class as \code{x} with reflected coordinates.
#'   
#' @param x an R object containing atomic coordinates.
#' @param p1 a numeric vector of length 3 containing the coordinates of the
#'   first point defining the reflexion plan. Can also be a 3x3 matrix or
#'   data.frame containing by row \code{p1}, \code{p2} and \code{p3}.
#' @param p2 a numeric vector of length 3 containing the coordinates of the
#'   second point defining the reflexion plane.
#' @param p3 a numeric vector of length 3 containing the coordinates of the
#'   third point defining the reflexion plane.
#' @param mask a logical vector indicating the set of coordinates to which to
#'   apply the reflexion.
#' @param cryst1 an object of class \sQuote{crystal} used to convert fractional
#'   into Cartesian coordinates (when needed).
#' @param \dots further arguments passed to or from other methods.
#'   
#' @seealso Helper functions for reflection with respect to a given Cartesian
#' plane or a plane defined by two lattice vectors:\cr \code{\link{Mxy}},
#' \code{\link{Myz}}, \code{\link{Mzx}}, \code{\link{Mab}}, \code{\link{Mbc}},
#' \code{\link{Mca}}\cr Passing from Cartesian to fractional coordinates (or Vis
#' Versa):\cr \code{\link{xyz2abc}}, \code{\link{abc2xyz}}
#' 
#' @examples 
#' # First lets read a pdb file
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' cell <- cell.coords(x)
#' visualize(x, mode = NULL)
#'
#' # Mirror operation with respect to the ab-plane
#' visualize(mirror(x, rep(0,3), p1=cell[, "a"], p2=cell[, "b"]), mode = NULL)
#' # Mirror operation with respect to the ab-plane for residue 1
#' visualize(mirror(x, rep(0,3), p1=cell[, "a"], p2=cell[, "b"],
#'    mask = x$atoms$resid == 1), mode = NULL)
#' 
#' @keywords manip
#'  
#' @name mirror
#' @export
mirror <- function(...)
  UseMethod("mirror")

#' @rdname mirror
#' @export
mirror.coords <- function(x, p1, p2 = NULL, p3 = NULL, mask = TRUE, cryst1 = NULL, ...){
  if(missing(p1))
    stop("Please specify at least 'p1'")
  if(is.null(p2) & is.null(p3)){
    if(ncol(p1)!=3 | nrow(p1)!=3)
      stop("When 'p2' and 'p3' are not specifyed, 'p1' must be a 3x3 matrix or data.frame")
    p3 <- p1[3,]
    p2 <- p1[2,]
    p1 <- p1[1,]
  } else {
    if(length(p3) != 3 | length(p2) != 3 | length(p1) != 3)
      stop("'p1', 'p2' and 'p3' must be vectors of length 3")
  }
  if(all(p1==p2)|all(p1==p3)|all(p2==p3))
    stop("'p1', 'p2' and 'p3' must be different to define the mirror")
  if(length(mask) != natom(x)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(x))
  }
  basis.ori <- basis(x)
  if(basis.ori != "xyz"){
    if(is.null(cryst1))
      stop("Please specify a 'crystal' obj to convert your fractional coordinates into Cartesian");
    x <- abc2xyz(x, cryst1 = cryst1)
  }
  v12 <- p2 - p1  
  v23 <- p3 - p2

  vn <- vectProd(v12,v23)
  vn <- vn/vectNorm(vn)
  x <- Txyz(x,  p1[1],  p1[2],  p1[3], mask=mask)
  rotM <- diag(3) - as.matrix(
    rbind(
      c(2*vn[1]*vn[1], 2*vn[1]*vn[2], 2*vn[1]*vn[3]),
      c(2*vn[2]*vn[1], 2*vn[2]*vn[2], 2*vn[2]*vn[3]),
      c(2*vn[3]*vn[1], 2*vn[3]*vn[2], 2*vn[3]*vn[3])
      )
    )
  x[mask,] <- coords(as.matrix(x[mask,])%*%rotM, basis = "xyz")
  x <- Txyz(x, -p1[1], -p1[2], -p1[3], mask=mask)
  return(x)
}

#' @rdname mirror
#' @export
mirror.pdb <- function(x, p1, p2 = NULL, p3 = NULL,
		mask = TRUE, cryst1 = x$crystal, ...) {
  coords(x) <- mirror(coords(x), p1=p1, p2=p2, p3=p3, mask=mask, cryst1=cryst1, ...)
  return(x)
}
