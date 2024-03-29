#' Create an Object of Class \sQuote{pdb}
#' 
#' Creates an object of class 'pdb'.
#' 
#' This function is the generic function to create objects of class
#' \sQuote{pdb}. The purpose of this class is to store the data of molecular
#' systems contained in PDB files. The default method of the \code{pdb} function
#' creates an object of class \sQuote{pdb} from its different components, i.e.:
#' \code{title}, \code{remark}, \code{cryst1}, \code{atoms} and \code{conect}. 
#' At least an object of class \sQuote{atoms} has to be specified. \cr\cr 
#' \code{is.pdb} tests if x is an object of class \sQuote{pdb}, i.e. if x has a
#' \dQuote{class} attribute equal to \code{pdb}.
#' 
#' @return 
#' \code{pdb} returns a list of class \sQuote{pdb} with the following components:
#' \item{title}{a character vector containing the TITLE records found in a PDB file.}
#' \item{remark}{a character vector containing the REMARK records found in a PDB file.}
#' \item{cryst1}{a list of class \sQuote{cryst1} containing the first CRYST1 record found in a PDB file. All others are ignored.}
#' \item{atoms}{a data.frame of class \sQuote{atoms} containing the ATOM and HETATM records found in a PDB file.}
#' \item{conect}{a data.frame of class \sQuote{conect} containing the CONECT records found in a PDB file.}
#' \cr
#' \code{is.pdb} returns TRUE if x is an object of class \sQuote{pdb} and FALSE otherwise.
#' 
#' @param atoms a data.frame of class \code{atoms} containing ATOM and HETATM records use to create the \code{pdb} object.
#' @param cryst1 a list of class \code{cryst1} containing the periodical bondary conditions and space group use to create the \code{pdb} object.
#' @param conect a data.frame of class \code{conect} containing CONECT records use to create the \code{pdb} object.
#' @param remark a character vector containing some REMARK records to be added to the \code{pdb} object.
#' @param title a character vector containing some TITLE records to be added to the \code{pdb} object.
#' @param x an R object to be tested
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso 
#' \code{\link{atoms}}, \code{\link{coords}}, \code{\link{cryst1}}, \code{\link{conect}} and \code{\link{read.pdb}}
#' 
#' @examples 
#' title  <- "This is just an example"
#' remark <- NULL
#' cryst1 <- cryst1(c(10,10,10))
#' atoms <- atoms(recname = c("ATOM","ATOM"), eleid = 1:2, elename = c("H","H"), alt = "",
#'                resname = c("H2","H2"), chainid = "", resid = c(1,1), insert = "",
#'                x1 = c(0,0), x2 = c(0,0), x3 = c(0,1), occ = c(0.0,0.0), temp = c(1.0,1.0),
#'                segid = c("H2","H2"))
#' conect <- conect(eleid.1 = c(1), eleid.2 = c(2))
#' x <- pdb(atoms = atoms, cryst1 = cryst1, conect = conect, remark = remark, title = title)
#' is.pdb(x)
#' 
#' @keywords classes
#' 
#' @name pdb
#' @export
pdb <- function(...)
  UseMethod("pdb")

#' @rdname pdb
#' @export
pdb.default <- function(atoms, cryst1 = NULL, conect = NULL, remark = NULL, title = NULL, ...)
{
  if(missing(atoms)) stop("Please specify at least an 'atoms' object")
  if(!is.atoms(atoms)) stop("'atoms' must be an object of class 'atoms'")
  
  if(!is.null(cryst1) & !is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  if(!is.null(conect) & !is.conect(conect)) stop("'conect' must be an object of class 'conect'")
  
  if(is.list(title ) | !is.null(dim(title ))) stop("'title' must be a vector of character strings")
  if(is.list(remark) | !is.null(dim(remark))) stop("'remark' must be a vector of character strings")
  
  if(!is.character(title ) & !is.null(title )) title  <- as.character(title )
  if(!is.character(remark) & !is.null(remark)) remark <- as.character(remark)
  
  to.return <- list(title = title, remark = remark, cryst1 = cryst1, atoms = atoms, conect = conect)
  
  class(to.return) <- c("pdb","list")
  return(to.return)
}

#' @rdname pdb
#' @export
is.pdb <- function(x)
{
  to.return <- any(class(x) == "pdb")
  return(to.return)
}
