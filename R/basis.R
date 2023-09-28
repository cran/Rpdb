#' The Basis of an Object
#' 
#' Functions to get or set the basis of an object containing atomic coordinates.
#' 
#' \code{basis} and \code{basis<-} are respectively generic accessor and 
#' replacement functions. The default methods get and set the \code{basis} 
#' attribute of an object containing atomic coordinates. This attribute indicate
#' the type basis vector used to express atomic coordinates.\cr \code{value} 
#' must be equal to \code{"xyz"}, for Cartesian, or \code{"abc"}, for fractional
#' coordinates.\cr The methods for objects of class \sQuote{pdb} get and set the
#' \code{basis} attribute of its \code{atoms} component.
#' 
#' @return \describe{ \item{For \code{basis}:}{NULL or a single element
#' character vector. (NULL is given if the object has no \code{basis}
#' attribute.)} \item{For \code{basis<-}:}{the updated object. (Note that the
#' value of \code{basis(x) <- value} is that of the assignment, value, not the
#' return value from the left-hand side.)} }
#' 
#' @param x an R object containing atomic coordinates.
#' @param value a single element character vector use to set the basis of 
#'   \code{x}.
#'
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}
#' 
#' @examples
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' basis(x)
#' x <- xyz2abc(x)
#' basis(x)
#' 
#' @keywords attribute
#'
#' @name basis
#' @export
basis <- function(x)
  UseMethod("basis")

#' @rdname basis
#' @export
basis.default <- function(x)
  attr(x, which = "basis")

#' @rdname basis
#' @export
'basis<-' <- function(x, value)
  UseMethod("basis<-", x)

#' @rdname basis
#' @export
'basis<-.default' <- function(x, value)
{
  if(!value %in% c("xyz","abc")) stop("Unrecognized basis")
  attr(x, "basis") <- value
  return(x)
}

#' @rdname basis
#' @export
basis.pdb <- function(x)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  to.return <- basis.default(x$atoms)
  return(to.return)
}

#' @rdname basis
#' @export
'basis<-.pdb' <- function(x, value)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  if(!value %in% c("xyz","abc")) stop("Unrecognized basis set")
  attr(x$atoms, "basis") <- value
  return(x)
}
