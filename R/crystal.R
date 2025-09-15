#' Create \sQuote{crystal} Object
#' 
#' Create an object of class \sQuote{crystal} containing the unit cell parameters
#' and the name of the space group to associate with an object of class 
#' \sQuote{pdb}. Note: the \sQuote{cryst1} class will be deprecated.
#' 
#' \code{crystal} is a generic function to create objects of class 
#' \sQuote{crystal}. The purpose of this class is to store CRYST1 records from 
#' PDB files which contain the unit cell parameters and the name of the space 
#' group of a molecular system stored in a PDB file. The default method of the 
#' \code{crystal} function creates an object of class \sQuote{crystal} from its 
#' different components, i.e.: \code{abc}, \code{abg} and \code{sgroup}. At 
#' least \code{abc} has to be specified. \cr\cr \code{is.crystal} tests if an 
#' object is of class \sQuote{crystal}, i.e. if it has a \dQuote{class} attribute 
#' equal to \code{crystal}.
#' 
#' @return Function \code{crystal} returns a list of class \sQuote{crystal} with the
#' following components:
#' \item{abc}{a numeric vector of length 3 containing the norms of the lattice
#'   vectors a, b and c.}
#' \item{abg}{a numeric vector of length 3 containing the angles between the
#'   lattice vectors \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.}
#' \item{sgroup}{a character string giving the Hermann-Mauguin symbol of the space group.}
#' 
#' Function \code{is.crystal} returns TRUE if \code{x} is an object of class \sQuote{crystal}
#'   and FALSE otherwise.
#' 
#' @param \dots further arguments passed to or from other methods.
#' @param abc a numeric vector of length 3 containing the norms of the lattice
#'   vectors a, b and c.
#' @param abg a numeric vector of length 3 containing the angles between the
#'   lattice vectors \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma}.
#' @param sgroup a character string giving the Hermann-Mauguin symbol of the
#'   space group.
#' @param x an R object to be tested.
#'   
#' @seealso  
#' \code{\link{cell.coords}}, \code{\link{pdb}}
#' 
#' @examples 
#' x <- crystal(abc = c(10, 10, 10), abg = c(90,90,90), sgroup = "P1")
#' is.crystal(x)
#'  
#' @keywords classes
#'  
#' @name crystal
#' @export
cryst1 <- function(...) {
	warning("Method will be deprecated!");
	UseMethod("crystal");
}
#' @export
crystal <- function(...)
	UseMethod("crystal");

# Helper: Deprecate cryst1;
checkArgCrystal = function(crystal, cryst1) {
	if( ! is.null(cryst1)) {
		if(is.null(crystal)) {
			warning("Arg cryst1 is deprecated!");
			return(cryst1);
		} else {
			stop("Please specify only crystal!",
				"Note: arg cryst1 will be deprecated!");
		}
	}
	return(crystal);
}

#' @rdname crystal
#' @export
crystal.default <- function(abc, abg = c(90, 90, 90), sgroup = "P1", ...)
{
	if(missing(abc)) stop("Please provide at leat 'abc'");
	if(is.crystal(abc)) {
		# TODO: What should be the default action?
		abc = abc$abc;
	}
	#
	to.return <- list(abc = abc, abg = abg, sgroup = sgroup)
	
	class(to.return) = c("crystal", "cryst1");
	return(to.return)
}

#' @rdname crystal
#' @export
is.cryst1 <- function(x)
{
  to.return <- inherits(x, c("crystal", "cryst1"));
  return(to.return)
}
#' @rdname crystal
#' @export
is.crystal <- function(x)
{
  to.return <- inherits(x, c("crystal", "cryst1"));
  return(to.return)
}


# Note: not yet exported;
as.crystal.character = function(x) {
	if(length(x) != 1) {
      warning("Multiple 'CRYST1' records have been found. Only the first record has been kept.")
      x <- x[1];
    }
    abc <- c(
      a = as.numeric(substr(x,  7, 15)),
      b = as.numeric(substr(x, 16, 24)),
      c = as.numeric(substr(x, 25, 33))
    )
    abg <- c(
      alpha = as.numeric(substr(x, 34, 40)),
      beta  = as.numeric(substr(x, 41, 47)),
      gamma = as.numeric(substr(x, 48, 54))
    )
    sgroup <- substr(x, 56, 66)
    
    if(any(is.na(abc))) warning("In 'crystal': 'abc' contains NA values");
    if(any(is.na(abg))) warning("In 'crystal': 'abg' contains NA values");
    if(sgroup == "") sgroup <- NULL

    crystal <- crystal.default(abc, abg, sgroup);
	return(crystal);
}
