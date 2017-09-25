#' Reinitialize Object Indexing
#' 
#' Reinitialize the indexing of an object.
#' 
#' \code{reindex} is a generic function to reinitialize the indexing of an object or its components.
#' The methods for objects of class \sQuote{atoms} reinitialize the residue and element IDs starting
#' from 1 and avoiding gaps in the indexes. For objects of class \sQuote{pdb} their \code{atoms} and
#' \code{conect} components are reindexed consistently.
#' 
#' @return Return an object of the same class as \code{x} with updated indexes.
#' 
#' @param x an R object.
#' @param eleid a single element logical vector indicating if elements IDs have to reindexed.
#' @param resid a single element logical vector indicating if residues IDs have to reindexed.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{pdb}}, \code{\link{atoms}}, \code{\link{subset.atoms}}, \code{\link{subset.pdb}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#' x <- subset(x, x$atoms$eleid \%in\% sample(x$atoms$eleid, 10))
#' print(x)
#' x <- reindex(x)
#' print(x)
#' 
#' @keywords manip
#' 
#' @name reindex
#' @export
reindex <- function(...)
  UseMethod("reindex")

#' @rdname reindex
#' @export
reindex.atoms <- function(x, eleid = TRUE, resid = TRUE, ...)
{
  if(eleid)
  {
    eleid <- as.factor(x$eleid)
    levels(eleid) <- 1:nlevels(eleid)
    x$eleid <- as.integer(as.character(eleid))
  }
  if(resid)
  {
    resid <- as.factor(x$resid)
    levels(resid) <- 1:nlevels(resid)
    x$resid <- as.integer(as.character(resid))
  }
  rownames(x) <- 1:nrow(x)
  return(x)
}

#' @rdname reindex
#' @export
reindex.pdb <- function(x, eleid = TRUE, resid = TRUE, ...)
{
  if(eleid)
  {
    if(anyDuplicated(x$atoms$eleid)){
      x$atoms$eleid <- 1:natom(x)
      if(!is.null(x$conect)){
        cat("Recalculating connectivity\n")
        x$conect <- conect(x)      
      }
    }
    else{
      if(!is.null(x$conect)){
        eleid.1 <- factor(x$conect$eleid.1, levels = x$atoms$eleid)
        eleid.2 <- factor(x$conect$eleid.2, levels = x$atoms$eleid)
        levels(eleid.1) <- 1:nlevels(eleid.1)
        levels(eleid.2) <- 1:nlevels(eleid.2)
        eleid.1 <- as.integer(as.character(eleid.1))
        eleid.2 <- as.integer(as.character(eleid.2))
        x$conect <- conect(eleid.1, eleid.2)
      }
      eleid <- factor(x$atoms$eleid, levels = x$atoms$eleid)
      neleid <- nlevels(eleid)
      levels(eleid) <- 1:neleid
      x$atoms$eleid <- as.integer(as.character(eleid))
    }
  }
  if(resid)
  {
    resid <- as.factor(x$atoms$resid)
    levels(resid) <- 1:nlevels(resid)
    x$atoms$resid <- as.integer(as.character(resid))
  }
  rownames(x$atoms) <- 1:nrow(x$atoms)
  return(x)
}
