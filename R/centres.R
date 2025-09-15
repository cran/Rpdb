#' Centres-of-Geometry and Centres-of-Mass
#' 
#' Computes centres-of-geometry and centres-of-mass of groups of atoms.
#' 
#' \code{centres} is a generic function to compute centres-of-geometry and 
#' centres-of-mass from an object containing atomic coordinates. For objects of 
#' class \sQuote{coords}, \sQuote{atoms} and \sQuote{pdb}, the coordinates of 
#' \code{x} are first splitted into groups defined by \code{factor} using the 
#' \code{\link{split}} function. For each group, the weighted mean of the 
#' \code{x1}, \code{x2} and \code{x3} components of \code{x} are calculated 
#' using \code{weights}. By default all atoms are assumed to have the same 
#' weight (calculation of centres-of-geometry). Finally, if \code{unplit = TRUE}
#' the coordinates of the centres are unsplitted using the \code{\link{unsplit}}
#' function to assign to each atom the coordinates of the centre to which they 
#' are attached (used for wrapping by groups). \cr \cr For objects of class 
#' \sQuote{atoms} and \sQuote{pdb} by default \code{factor} is set to 
#' \code{x$resid} and \code{x$coordinates$resid}, respectively, to compute the 
#' centre-of-geometry of the different resdiues. Notice that coordinates can be 
#' neglected for the calculation of the centres using NA values in 
#' \code{factor}.
#' 
#' @return Return an object of class \sQuote{coords} containing the coordinates
#'   of centres.
#'   
#' @param x an R object containing atomic coordinates.
#' @param factor a factor used to split the atomic coordinates by groups to
#'   compute multiple centres.
#' @param weights a numerical vector containing atomic weights used to compute
#'   centres-of-mass.
#' @param unsplit a logical value indicating whether the coordinates of the
#'   centres have to be unsplit to repeat their coordinates for each atom used
#'   for their calculation (used for wrapping by groups).
#' @param na.rm a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}},
#' \code{\link{elements}} \cr\cr and \code{\link{split}}, \code{\link{unsplit}},
#' \code{\link[base]{factor}} for details about splitting data sets.
#' 
#' @examples
#' # First lets read a pdb file
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' 
#' # Centres-of-geometry of the residues
#' centres(x)
#' 
#' # Centre-of-geometry of the whole structure
#' centres(x, factor = rep(1, natom(x)))
#' # or
#' centres(coords(x))
#' 
#' # Centres-of-geometry of the PCB and DCB residues
#' centres(x, factor = x$atoms$resname)
#' 
#' # Knowing the name of the elements forming
#' # the C60 of the PCBM molecules (PCB residues)
#' # we can compute the centres-of-geometry of
#' # the C60 by neglecting the other atoms of the
#' # PCB residues.
#' C60.elename <- paste0("C", sprintf("%0.3d", 1:60))
#' 
#' is.PCB <- x$atoms$resname == "PCB" # Produce a mask to select only the PCB residues
#' is.C60 <- is.PCB & x$atoms$elename %in% C60.elename # Produce a mask to keep only the C60
#' 
#' F <- x$atoms$resid # We use the residue IDs to split the coordinates
#' F[!is.C60] <- NA # We keep only the atoms of the C60
#' 
#' C60.centres <- centres(x, factor = F)
#' 
#' # Lets check the position of the C60 centres
#' visualize(x , mode = NULL)
#' spheres3d(C60.centres)
#' text3d(Ty(C60.centres, 2), text=paste0("PCB_", rownames(C60.centres)), cex=2)
#' 
#' # Centres-of-mass of the resdiues
#' symb <- toSymbols(x$atoms$elename) # Convert elename into elemental symbols
#' # Find the mass of the element in the periodic table
#' w <- elements[match(symb, elements[,"symb"]), "mass"] 
#' centres(x, weights = w)
#' 
#' @keywords manip 
#' 
#' @name centres
#' @export
centres <- function(x, ...)
  UseMethod("centres")

#' @rdname centres
#' @export
centres.coords <- function(x, factor = NULL, weights = NULL, unsplit = FALSE, na.rm = FALSE, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")

  if(is.null(factor )) factor  <- rep("", natom(x))
  if(is.null(weights)) weights <- rep(1 , natom(x))

  w  <- split(weights, factor)
  x1 <- split(x$x1   , factor)
  x2 <- split(x$x2   , factor)
  x3 <- split(x$x3   , factor)
  
  w.mean <- function(x, w)
    sum(x*w/sum(w, na.rm = na.rm))
  
  x1.mean <- mapply(w.mean, x1, w, SIMPLIFY = FALSE)
  x2.mean <- mapply(w.mean, x2, w, SIMPLIFY = FALSE)
  x3.mean <- mapply(w.mean, x3, w, SIMPLIFY = FALSE)
  
  if(unsplit){
    x1.mean <- unsplit(x1.mean, factor)
    x2.mean <- unsplit(x2.mean, factor)
    x3.mean <- unsplit(x3.mean, factor)
  }
  
  x1.mean <- unlist(x1.mean)
  x2.mean <- unlist(x2.mean)
  x3.mean <- unlist(x3.mean)
  
  to.return <- coords.default(x1.mean, x2.mean, x3.mean, basis(x))
  
  return(to.return)
}

#' @rdname centres
#' @export
centres.atoms <- function(x, factor = NULL, weights = NULL, unsplit = FALSE, na.rm = FALSE, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  if(is.null(factor)) factor <- x$resid
  
  to.return <- centres.coords(coords(x), factor, weights, unsplit, na.rm)
  return(to.return)
}

#' @rdname centres
#' @export
centres.pdb <- function(x, factor = NULL, weights = NULL, unsplit = FALSE, na.rm = FALSE, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  to.return <- centres.atoms(x$atoms, factor, weights, unsplit, na.rm)
  return(to.return)
}
