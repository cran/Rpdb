#' Merging Molecular Systems
#' 
#' Merge two objects contaning atomic coordinates
#' 
#' To merge \code{x} and \code{y} they must have the same \code{basis} 
#' attributes (see \code{\link{basis}}). \cr\cr For objects of class 
#' \sQuote{coords} and \sQuote{atoms} the atomic coordinates are directly merged
#' by row. \cr\cr For objects of class \sQuote{pdb}, the \code{atoms} and 
#' \code{conect} components of the two \code{pdb} objects are merged by row and 
#' the \code{crystal} component of \code{x} is used to build the returned 
#' object. \cr\cr For objects of class \sQuote{atoms} and \sQuote{pdb} the 
#' residue and element IDs of \code{y} are shifted to avoid any confusion with 
#' those of \code{x}. If \code{reindex == TRUE} the \code{\link{reindex}} function
#' is called to reinitialize the indexing of the returned object.
#' 
#' @return Return an object of the same class as \code{x} and \code{y} merging 
#'   \code{x} and \code{y}. If \code{x} and \code{y} have different \code{basis}
#'   attributes an error is returned.
#'   
#' @param x,y objects of class \code{coords} to be merged.
#' @param reindex a single element logical vector indicating if residue and
#'   element IDs have to be reindexed after merging.
#' @param \dots further arguments passed to or from other methods.
#'   
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}},
#'   \code{\link{basis}}, \code{merge}, \code{merge.data.frame}
#'   
#' @examples 
#' c1 <- coords( 1:3 ,  4:6 ,  7:9 , basis = "xyz")
#' c2 <- coords(10:12, 13:15, 16:18, basis = "xyz")
#' merge(c1,c2)
#' 
#' \donttest{
#' ## Merging objects with different basis sets returns an error.
#' c2 <- coords(9:11, 12:14, 15:17, basis = "abc")
#' try(merge(c1,c2))
#' }
#' 
#' ## Prepare a Pentacene/C70 dimer
#' C70 <- read.pdb(system.file("examples/C70.pdb", package="Rpdb"))
#' Pen <- read.pdb(system.file("examples/Pentacene.pdb", package="Rpdb"))
#' x <- merge(Tz(C70, 3.5, thickness=0.5), Pen)
#'   
#' @keywords manip
#' 
#' @name merge.coords
#' @export
merge.coords <- function(x, y, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.coords(y)) stop("'y' must be an object of class 'coords'")
  
  if(basis(x) != basis(y)) stop("'x' and 'y' basis differ")
  
  to.return <- rbind(x,y)
  return(to.return)
}

#' @rdname merge.coords
#' @export
merge.atoms <- function(x, y, reindex = TRUE, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  if(!is.atoms(y)) stop("'y' must be an object of class 'atoms'")
  
  if(basis(x) != basis(y)) stop("'x' and 'y' basis differ")
  
  y.resid <- as.integer(y$resid + max(x$resid))
  y.eleid <- as.integer(y$eleid + max(x$eleid))
  
  y$resid <- y.resid
  y$eleid <- y.eleid
  
  to.return <- rbind(x,y)
  if(reindex) to.return <- reindex.atoms(to.return)
  
  return(to.return)
}

#' @rdname merge.coords
#' @export
merge.pdb <- function(x, y, reindex = TRUE, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  if(!is.pdb(y)) stop("'y' must be an object of class 'pdb'")
  
  if(basis(x) != basis(y)) stop("'x' and 'y' basis differ")
  
	xcryst = x$crystal; ycryst = y$crystal;
	if( any(xcryst$abc != ycryst$abc) ||
		any(xcryst$abg != ycryst$abg) ||
		any(xcryst$sgroup != ycryst$sgroup))
		warning("Different periodical boundary conditions are defined for 'x' and 'y'.",
			"'x$crystal' has been kept.");
  
	title  = unique(c(x$title , y$title ));
	remark = unique(c(x$remark, y$remark));
	cryst1 = x$crystal;
  
  y$conect$eleid.1 <- y$conect$eleid.1 + max(x$atoms$eleid)
  y$conect$eleid.2 <- y$conect$eleid.2 + max(x$atoms$eleid)

  eleid.1 <- c(x$conect$eleid.1, y$conect$eleid.1)
  eleid.2 <- c(x$conect$eleid.2, y$conect$eleid.2)
  conect      <- conect(eleid.1, eleid.2)

  atoms <- merge.atoms(x$atoms, y$atoms, reindex = FALSE)
  
  to.return   <- pdb(atoms, cryst1, conect, remark, title)
  if(reindex) to.return <- reindex.pdb(to.return)
  
  return(to.return)
}

