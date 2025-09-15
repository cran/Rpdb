#' Create \sQuote{atoms} Object
#' 
#' Creates an object of class \sQuote{atoms} containing the data related to ATOM
#' and HETATM records of a PDB file.
#' 
#' \code{atoms} is a generic function to create objects of class \sQuote{atoms}.
#' The purpose of this class is to store ATOM and HETATM records from PDB files.
#' The default method creates a \code{atoms} object from its different
#' components, i.e.: \code{recname}, \code{eleid}, \code{elename}, \code{alt},
#' \code{resname}, \code{chainid}, \code{resid}, \code{insert}, \code{x1},
#' \code{x2}, \code{x3}, \code{occ}, \code{temp}, \code{segid} and \code{basis}.
#' All the arguments have to be specified except \code{basis} which by default
#' is set to "xyz" (Cartesian coordinates).\cr \code{is.atoms} tests if an object
#' of class \sQuote{atoms}, i.e. if it has a \dQuote{class} attribute equal to
#' \code{atoms}.
#' 
#' @return 
#' \code{atoms} returns a data.frame of class \sQuote{atoms} with the following components:
#' \describe{
#' \item{recname}{a character vector containing the record name for each element.}
#' \item{eleid}{a integer vector containing the element ID for each element.}
#' \item{elename}{a character vector containing the element name for each element.}
#' \item{alt}{a character vector containing the alternate location indicator for each element.}
#' \item{resname}{a character vector containing the residue name for each element.}
#' \item{chainid}{a character vector containing the chain ID for each element.}
#' \item{resid}{a integer vector containing the residue ID for each element.}
#' \item{insert}{a character vector containing the codes for insertion of residue for each element.}
#' \item{x1, x2, x3}{a numeric vector containing the first, second and third coordinate for each element.}
#' \item{occ}{a numeric vector containing the occupencie for each element.}
#' \item{temp}{a numeric vector containing the temperature factor for each element.}
#' \item{segid}{a character vector containing the segment ID for each element.}
#' \item{basis}{a single element character vector indicating the type of basis vector used to express the atomic coordinates.}
#' }
#' \code{is.atoms} returns TRUE if \code{x} is an object of class \sQuote{atoms} and FALSE otherwise.
#' 
#' @param \dots arguments passed to methods.
#' @param recname a character vector containing the record name for each element.
#' @param eleid a integer vector containing the element ID for each element.
#' @param elename a character vector containing the element name for each element.
#' @param alt a character vector containing the alternate location indicator for each element.
#' @param resname a character vector containing the residue name for each element.
#' @param chainid a character vector containing the chain ID for each element.
#' @param resid a integer vector containing the residue ID for each element.
#' @param insert a character vector containing the codes for insertion of residue of each element.
#' @param x1,x2,x3 a numeric vector containing the first, second and third coordinate for each element.
#' @param occ a numeric vector containing the occupancie for each element.
#' @param temp a numeric vector containing the temperature factor for each element.
#' @param segid a character vector containing the segment ID for each element.
#' @param basis a single element character vector indicating the type of basis vector used to express the atomic coordinates.
#' @param x an R obecjt to be tested.
#' 
#' @seealso \code{\link{basis}}, \code{\link{coords}}, \code{\link{pdb}}
#' 
#' @examples 
#' x <- atoms(recname = c("ATOM","ATOM"), eleid = 1:2, elename = c("H","H"), alt = "",
#'   resname = c("H2","H2"), chainid = "", resid = c(1,1), insert = "",
#'   x1 = c(0,0), x2 = c(0,0), x3 = c(0,1), occ = c(0.0,0.0), temp = c(1.0,1.0),
#'   segid = c("H2","H2"))
#' print(x)
#' is.atoms(x)
#' 
#' @keywords classes
#' 
#' @import rgl
#' @importFrom grDevices rgb
#' 
#' @name atoms
#' @export
atoms <- function(...)
  UseMethod("atoms")

#' @rdname atoms
#' @export
atoms.default <- function(recname, eleid, elename, alt,
                                resname, chainid, resid, insert,
                                x1, x2, x3, occ, temp, segid, basis = "xyz", ...)
{
  
  recname <- as.character(recname)
  eleid   <- suppressWarnings(as.integer(eleid))
  elename <- as.character(elename)
  alt     <- as.character(alt)
  resname <- as.character(resname)
  chainid <- as.character(chainid)
  resid   <- suppressWarnings(as.integer(resid))
  insert  <- as.character(insert)
  x1      <- suppressWarnings(as.numeric(x1))
  x2      <- suppressWarnings(as.numeric(x2))
  x3      <- suppressWarnings(as.numeric(x3))
  occ     <- suppressWarnings(as.numeric(occ))
  temp    <- suppressWarnings(as.numeric(temp))
  segid   <- as.character(segid)
  
  if(any(is.na(eleid))) warning("In 'atoms': 'eleid' contains NA values")
  if(any(is.na(resid))) warning("In 'atoms': 'resid' contains NA values")
  if(any(is.na(   x1))) warning("In 'atoms':    'x1' contains NA values")
  if(any(is.na(   x2))) warning("In 'atoms':    'x2' contains NA values")
  if(any(is.na(   x3))) warning("In 'atoms':    'x3' contains NA values")
  if(any(is.na(  occ))) warning("In 'atoms':   'occ' contains NA values")
  if(any(is.na( temp))) warning("In 'atoms':  'temp' contains NA values")
  
  atoms <- data.frame(
    recname = recname,
    eleid   = eleid  ,
    elename = elename,
    alt     = alt    ,
    resname = resname,
    chainid = chainid,
    resid   = resid  ,
    insert  = insert ,
    x1     = x1      ,
    x2     = x2      ,
    x3     = x3      ,
    occ     = occ    ,
    temp    = temp   ,
    segid   = segid  ,
    stringsAsFactors = FALSE
  )
  attr(atoms, "basis") <- basis
  
  class(atoms) <- c("atoms","coords","data.frame")
  return(atoms)
}

#' @rdname atoms
#' @export
is.atoms <- function(x)
{
  to.return <- any(attr(x,which="class") == "atoms")
  return(to.return)
}
