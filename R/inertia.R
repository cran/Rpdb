#' Moment of Inertia of a Molecular System
#' 
#' Computes the inertia tensor of a molecular system from atomic coordinates and
#' masses.
#' 
#' \code{inertia} is a generic function to compute the inertia tensor of a
#' molecular system. For object of class \sQuote{coords} both atomic coordinates
#' and masses have to be speifyed. For object of class \sQuote{atoms} the masses
#' are determined from the \code{elename} component of the object (see
#' \code{\link{toSymbols}} and \code{\link{masses}}). For object of class
#' \sQuote{pdb} the \code{atoms} component is used.
#' 
#' @return Return the inertia tensor in a 3x3 matrix.
#' 
#' @param x an R object containing atomic coordinates.
#' @param m a numeric vector containing atomic masses.
#' @param \dots further arguments passed to or from other methods.
#'
#' @seealso 
#' \code{\link{toSymbols}}, \code{\link{masses}}, \code{\link{viewInertia}}
#'
#' @examples 
#' C70 <- read.pdb(system.file("examples/C70.pdb",package="Rpdb"))
#' inertia(C70)
#' visualize(C70, mode = NULL)
#' viewXY()
#' viewInertia(C70)
#' 
#' @keywords manip
#' 
#' @name inertia
#' @export
inertia <- function(...)
  UseMethod("inertia")

#' @rdname inertia
#' @export
inertia.coords <- function(x, m = NULL, ...){
  if(!is.coords(x))
    stop("'x' must be an object of class coords")
  if(is.null(m))
    stop("Please specify the masses")
  if(any(is.na(m))|any(is.na(x)))
    stop("NA values not permetted")
  Ixx<-sum(m*(x$x2^2+x$x3^2))
  Iyy<-sum(m*(x$x1^2+x$x3^2))
  Izz<-sum(m*(x$x1^2+x$x2^2))
  Ixy<-sum(m*(x$x1*x$x2))
  Ixz<-sum(m*(x$x1*x$x3))
  Iyz<-sum(m*(x$x2*x$x3))
  I<-matrix(c(Ixx,-Ixy,-Ixz,-Ixy,Iyy,-Iyz,-Ixz,-Iyz,Izz),ncol=3)
  return(I)
}

#' @rdname inertia
#' @export
inertia.atoms <- function(x, m = NULL, ...){
  if(is.null(m))
    m <- masses(toSymbols(x$elename))
  inertia(coords(x), m)
}

#' @rdname inertia
#' @export
inertia.pdb <- function(x, m = NULL, ...)
  inertia(x$atoms, m)
