#' The Atomic Coordinates of an Object
#' 
#' Get or set the atomic coordinates (either Cartesian or fractional 
#' coordinates) of an object.
#' 
#' The purpose of the \sQuote{coords} class is to store the coordinates of a 
#' molecular system and facilitate their manipulation when passing from the 
#' Cartesian to fractional coordinates and vice versa.\cr
#' \code{coords} and \code{coords<-} are generic accessor and replacement
#' functions. \code{as.coords} is an alias to \code{coords}.\cr
#' The default method of the \code{coords} function is actually a builder allowing 
#' to create a \sQuote{coords} object from its different components, i.e.: 
#' \code{x1}, \code{x2}, \code{x3}, and \code{basis}. All the arguments have to 
#' be specified except 'basis' which by default is set to "xyz" (Cartesian 
#' coordinates). \cr\cr
#' For an object of class \sQuote{atoms}, the accessor function
#' extracts its \code{x1}, \code{x2} and \code{x3} components as well 
#' as its \code{basis} attribute to create a \sQuote{coords} object.
#'
#' The replacement function sets \code{x1}, \code{x2} and \code{x3} components,
#' as well as the \code{basis} attribute. \cr\cr
#' For an object of class \sQuote{coords}, the accessor function
#' returns the \sQuote{coords} object as is.
#'
#' For an object of class \sQuote{pdb}, the accessor function extracts the \code{x1},
#' \code{x2} and \code{x3} components, as well as the \code{basis} attribute of its 
#' \code{atoms} component to create a \sQuote{coords} object. The replacement 
#' function sets the \code{x1}, \code{x2} and \code{x3} components as well as the
#' \code{basis} attribute of its \code{atoms} component. \cr\cr
#'
#' For \sQuote{matrix} and \sQuote{data.frame} objects, when \code{basis == NULL}
#' this function searches x, y, z or a, b, c columns in \code{x}.\cr
#' If x, y, z columns are found, they are used to set the first, second and
#' third coordinates of the returned \sQuote{coords} object.
#' In that case the basis set of \code{x} is set to \code{"xyz"}.\cr
#' If a, b, c columns are found they are used to a 
#' set the first, second and third coordinates of the returned \sQuote{coords} 
#' object. In that case the basis set of \code{x} is set to \code{"abc"}.\cr If 
#' the function doesn't found neither the x, y, z nor the a, b, c columns an 
#' error is returned.\cr When \code{basis!=NULL} it has to be equal to 
#' \code{"xyz"} or \code{"abc"} and \code{x} must have exactly 3 columns. \cr\cr
#' \code{is.coords} tests if x is an object of class \sQuote{coords}, i.e. if x 
#' has a \dQuote{class} attribute equal to \code{coords}.
#' 
#' @return The accessor function returns a data.frame of class \sQuote{coords}
#' whose columns contain the three coordinates of the atoms of a molecular
#' system. The coordinates can either be Cartesian (\code{basis} attribute equal
#' to \code{"xyz"}) or fractional coordinates (\code{basis} attribute equal to
#' \code{"abc"}). \cr\cr The replacement function returns an object of the same
#' class as \code{x} with updated coordinates. \cr\cr \code{is.coords} returns
#' TRUE if x is an object of class \sQuote{coords} and FALSE otherwise
#' 
#' @param x1,x2,x3 numeric vectors containing the first, second and third coordinates.
#' @param basis a single element character vector indicating the type of basis vector
#'   used to express the atomic coordinates.
#' @param x an R object containing atomic coordinates.
#' @param value an object of class \sQuote{coords} used for replacement
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso 
#' \code{\link{basis}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' is.coords(x)
#' is.coords(x$atoms)
#' 
#' ## Replace the coordinates of x by translated coordinates
#' coords(x) <- coords(Tz(x, 10))
#' coords(x)
#' 
#' @keywords classes manip
#' 
#' @name coords
#' @export
coords <- function(...)
  UseMethod("coords")

#' @rdname coords
#' @export
'coords<-' = function(x, value)
  UseMethod("coords<-", x)

#' @rdname coords
#' @export
'as.coords' = function(...) {
	UseMethod("coords");
}

#' @rdname coords
#' @export
coords.default <- function(x1, x2, x3, basis = "xyz", ...)
{
  if(!basis %in% c("xyz", "abc")) stop("Unrecognized basis")
  
  to.return <- data.frame(x1,x2,x3)
  attr(to.return, which = "basis") <- basis
  
  class(to.return) <- c("coords", "data.frame")
  
  return(to.return)
}

#' @rdname coords
#' @export
coords.data.frame <- function(x, basis = NULL, ...)
{
  if(!is.data.frame(x)) stop("'x' must be a 'data.frame'")
  
  if(is.null(basis)){
    if(all(c("x","y","z") %in% names(x))){
      x <- x[, c("x","y","z")];
      basis <- "xyz";
    }
    else if(all(c("a","b","c") %in% names(x))){
      x <- x[,c("a","b","c")]
      basis <- "abc"
    }
    else stop("Can not convert this 'data.frame' into 'coords': Coordinates not found")
  }
  else if( ! basis %in% c("xyz","abc")) stop("Unrecognized 'basis'");
  if(ncol(x) != 3L) stop("'x' must be a three-columns data.frame")

  to.return <- coords.default(x[,1], x[,2], x[,3], basis = basis)

  return(to.return)
}

#' @rdname coords
#' @export
coords.matrix <- function(x, basis = NULL, ...){
  if( ! is.matrix(x)) stop("'x' must be a 'matrix'");
  
  to.return <- coords.data.frame(as.data.frame(x), basis = basis, ...)
  return(to.return)
}

#' @rdname coords
#' @export
coords.coords <- function(x, ...) {
  if( ! is.coords(x)) stop("'x' must be a 'coords'-object");
  return(x);
}

#' @rdname coords
#' @export
coords.atoms <- function(x, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  to.return <- coords(x$x1,x$x2,x$x3,basis.default(x))
  return(to.return)
}

#' @rdname coords
#' @export
'coords<-.atoms' <- function(x, value)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  if(!is.coords(value)) stop("'value' must be an object of class 'coords'")
  if(nrow(x) != nrow(value)) stop(paste("arguments imply different number of rows: ",nrow(x),", ",nrow(value),sep=""))
  x[c("x1","x2","x3")] <- value
  basis(x) <- basis.default(value)
  
  return(x)
}

#' @rdname coords
#' @export
coords.pdb <- function(x, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  to.return <- coords.atoms(x$atoms)
  return(to.return)
}

#' @rdname coords
#' @export
'coords<-.pdb' <- function(x, value)
{
  if( ! is.pdb(x)) stop("'x' must be an object of class 'pdb'");
  
  if( ! is.coords(value)) stop("'value' must be an object of class 'coords'");
  if(nrow(x$atoms) != nrow(value)) {
		stop(paste("arguments imply different number of rows: ",
			nrow(x$atoms), ", ", nrow(value), sep=""));
	}
  coords(x$atoms) <- value
  
  return(x)
}

#' @rdname coords
#' @export
is.coords <- function(x)
{
  to.return <- inherits(x, "coords");
  return(to.return)
}
