#' Mass of Chemical Elements
#' 
#' Determine the mass of chemical elements
#' 
#' \code{masses} is a generic function to determine the mass of chemical 
#' elements. \cr\cr For objects of class \sQuote{pdb}: \itemize{ \item First the
#' element names are converted into element symbols using the \code{toSymbols} 
#' function. \item Then their masses are taken from the \code{elements} data 
#' set. } \code{NA} values are returned for unrecognized elements.
#' 
#' @return 
#' Return a numeric vector containing the mass of chemical elements.
#' 
#' @param \dots further arguments passed to or from other methods.
#' @param x either a character or an integer vector containing element symbols
#'   or atomic numbers, or an object of class \sQuote{pdb} from which element
#'   symbols are determined (see details).
#'   
#' @seealso 
#' \code{\link{toSymbols}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#' masses(x)
#' 
#' masses(c("C","Cl",NA,"AA","N"))
#' 
#' @keywords manip
#'       
#' @name masses
#' @export
masses <- function(...)
  UseMethod("masses")

#' @rdname masses
#' @export
masses.default <- function(x, ...) {
  if(is.character(x)) return(Rpdb::elements$mass[match(x, Rpdb::elements$symb)])
  else if(is.numeric(x) & x == round(x) ) return(Rpdb::elements$mass[match(x, Rpdb::elements$num)])
  else stop("Bad argument: 'x' must be a character or an integer vector")
}

#' @rdname masses
#' @export
masses.pdb <- function(x, ...) {
  x <- toSymbols(x$atoms$elename)
  return(masses(x))
}