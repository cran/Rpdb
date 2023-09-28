#' Properties of a Unit Cell
#' 
#' Compute the Cartesian coordinates of lattice vectors, the volume or the density of a unit cell.
#'
#' \code{cell.coords} is a generic function which computes a 3x3 matrix whose columns contrain the Cartesian coordinates of lattice vectors.
#' The 'a' and 'b' vectors are assumed to be respectively along the x-axis and in the xy-plane.
#' The default method takes directly the lattice parameters as arguments.
#' For objects of class \code{\link{cryst1}} the lattice parameters are first extracted from the object and then the default method is called.
#' For objects of class \code{\link{pdb}} the lattice parameters are extracted from their \code{cryst1} component and the default method is called.
#' \cr
#' \cr
#' \code{cell.volume} is a generic function to compute the volume of a unit cell.
#' For objects of class \sQuote{cryst1}, the unit cell parameters are directly used to compute the volume.
#' For objects of class \sQuote{pdb}, their \code{cryst1} component is used.
#' \cr
#' \cr
#' \code{cell.density} is a generic function to compute the density of a unit cell.
#' For objects of class \sQuote{pdb}:
#' First the volume of the unit cell is calculated by calling the \code{cell.volume} function on the \code{cryst1} component of the \sQuote{pdb} object.
#' Then the element names are converted into element symbols using the \code{toSymbols} function and their masses are taken from the \code{elements} data set.
#' Finally the density is calculated using the sum of the atomic masses and the volume of the unit cell.
#'
#' @return 
#' \code{cell.coords} returns a 3x3 matrix containing the Cartesian coordinates of lattice vectors arranged by columns.\cr
#' \code{cell.volume} returns a single element numeric vector containing the volume of the unit cell in Angstrom cube.\cr
#' \code{cell.density} returns a single element numeric vector containing the density of the unit cell in g.cm-3.
#'
#' @param abc a length 3 numeric vector containing the length of the a, b and c lattice vectors.
#' @param abg a length 3 numeric vector containing the angles (degrees) between the a, b and c lattice vectors (alpha, beta, gamma).
#' @param digits an integer used to round the lattice vectors coordinates.
#' @param x an R object containing lattice parameters.
#' @param masses a numeric vector containing atomic masses.
#' @param volume a single element numeric vector containing the volume of the unit cell in Angstrom cube.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{cryst1}}, \code{\link{pdb}}, \code{\link{xyz2abc}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' cell.volume(x)
#' cell.density(x)
#' cell.coords(x)
#'
#' @keywords manip
#'
#' @name cellProperties
#' @export
cell.coords <- function(...)
  UseMethod("cell.coords")

#' @rdname cellProperties
#' @export
cell.coords.default <- function(abc, abg = c(90,90,90), digits = 3, ...)
{
  if(missing(abc)) stop("Please provide at list a 'abc' vector containing the length of the lattice vectors")
  if(length(abc) != 3) stop("'abc' must be a vector of length 3")
  if(length(abg) != 3) stop("'abg' must be a vector of length 3")
  
  abg <- abg*pi/180
  
  M <- matrix(ncol=3,nrow=3)
  M[ ,1] <- c(abc[1],0,0)
  M[ ,2] <- c(abc[2]*cos(abg[3]),abc[2]*sin(abg[3]),0)
  M[1,3] <-   abc[3]*(cos(abg[2]))
  M[2,3] <-   abc[3]*(cos(abg[1])-cos(abg[2])*cos(abg[3]))/sin(abg[3])
  M[3,3] <-   abc[3]*sqrt(1+2*cos(abg[1])*cos(abg[2])*cos(abg[3])-(cos(abg[1]))^2-(cos(abg[2]))^2-(cos(abg[3]))^2)/sin(abg[3])
  M <- round(M, digits=3)
  dimnames(M) <- list(c("x","y","z"), c("a","b","c"))
  
  return(M)
}

#' @rdname cellProperties
#' @export
cell.coords.cryst1 <- function(x, digits = 3, ...)
{
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
  
  M <- cell.coords.default(x$abc, x$abg, digits)
  return(M)
}

#' @rdname cellProperties
#' @export
cell.coords.pdb <- function(x, digits = 3, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'atoms'")
  if(is.null(x$cryst1)) stop("'x' must contained a 'cryst1' object")
  
  M <- cell.coords.cryst1(x$cryst1, digits = 3)
  return(M)
}

#' @rdname cellProperties
#' @export
cell.volume <- function(...)
  UseMethod("cell.volume")

#' @rdname cellProperties
#' @export
cell.volume.cryst1 <- function(x, ...)
{
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
  
  V <-  prod(x$abc)*sqrt(1 - sum(cos(x$abg*pi/180)^2) + 2*prod(cos(x$abg*pi/180)))
  #   attr(V, "unit") <- "AngtromCube"
  return(V)
}

#' @rdname cellProperties
#' @export
cell.volume.pdb <- function(x, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'atoms'")
  if(is.null(x$cryst1)) stop("'x' must contained a 'cryst1' object")
  
  V <- cell.volume.cryst1(x$cryst1)
  return(V)
}

#' @rdname cellProperties
#' @export
cell.density <- function(...)
  UseMethod("cell.density")

#' @rdname cellProperties
#' @export
cell.density.default <- function(masses, volume, ...) {
  Na <- Rpdb::universalConstants["Na","Value"]
  d <- sum(masses)/(volume*1E-24*Na)
  #   attr(d, "unit") <- "g.cm-3"
  return(d)
}

#' @rdname cellProperties
#' @export
cell.density.pdb <- function(x, ...) {
  M <- masses(x)
  V <- cell.volume(x)
  d <- cell.density(M, V)
  return(d)
}
