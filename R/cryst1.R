#' Create \sQuote{cryst1} Object
#' 
#' Create an object of class \sQuote{cryst1} containing the unit cell parameters
#' and the name of the space group to associate with an object of class 
#' \sQuote{pdb}.
#' 
#' \code{cryst1} is a generic function to create objects of class 
#' \sQuote{cryst1}. The purpose of this class is to store CRYST1 records from 
#' PDB files which contain the unit cell parameters and the name of the space 
#' group of a molecular system stored in a PDB file. The default method of the 
#' \code{cryst1} function create an object of class \sQuote{cryst1} from its 
#' different components, i.e.: \code{abc}, \code{abg} and \code{sgroup}. At 
#' least \code{abc} has to be specified. \cr\cr \code{is.cryst1} tests if an 
#' object of class \sQuote{cryst1}, i.e. if it has a \dQuote{class} attribute 
#' equal to \code{cryst1}.
#' 
#' @return \code{cryst1} returns a list of class \sQuote{cryst1} with the
#' following components: \item{abc}{a numeric vector of length 3 containing the
#' norms of the lattice vectors a, b and c.} \item{abg}{a numeric vector of
#' length 3 containing the angles between the lattice vectors \eqn{\alpha},
#' \eqn{\beta} and \eqn{\gamma}.} \item{sgroup}{a character string giving the
#' Hermann-Mauguin symbol of the space group.} \cr\cr \code{is.Cryst1} returns
#' TRUE if \code{x} is an object of class \sQuote{cryst1} and FALSE otherwise.
#' 
#' @param \dots further arguments passed to or from other methods.
#' @param abc a numeric vector of length 3 containing the norms of the lattice
#'   vectors a, b and c.
#' @param abg a numeric vector of length 3 containing the angles between the
#'   lattice vectors \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.
#' @param sgroup a character string giving the Hermann-Mauguin symbol of the
#'   space group.
#' @param x an R obecjt to be tested.
#'   
#' @seealso  
#' \code{\link{cell.coords}}, \code{\link{pdb}}
#' 
#' @examples 
#' x <- cryst1(abc = c(10, 10, 10), abg = c(90,90,90), sgroup = "P1")
#' is.cryst1(x)
#'  
#' @keywords classes
#'  
#' @name cryst1
#' @export
cryst1 <- function(...)
  UseMethod("cryst1")

#' @rdname cryst1
#' @export
cryst1.default <- function(abc, abg = c(90, 90, 90), sgroup = "P1", ...)
{
  if(missing(abc)) stop("Please provide at leat 'abc'")
  to.return <- list(abc = abc, abg = abg, sgroup = sgroup)
  
  class(to.return) <- "cryst1"
  return(to.return)
}

#' @rdname cryst1
#' @export
is.cryst1 <- function(x)
{
  to.return <- any(attr(x,which="class") == "cryst1")
  return(to.return)
}