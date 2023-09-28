#' Subsetting \sQuote{atoms} and \sQuote{pdb} Objects
#' 
#' Return subsets of \sQuote{atoms} or \sQuote{pdb} objects which meet
#' conditions.
#' 
#' For an \sQuote{atoms} object the method is similar to the data.frame method
#' (see \code{\link{subset}}) but allow to directly reindex the elements and
#' residues IDs. For a \sQuote{pdb} object subsetting is applied on the
#' \code{atoms} and \code{conect} components of the object in a consistent way.
#' First the \code{atoms} component is subsetted and then the \code{conect}
#' component is filtered to keep only the conectivity for the subset.
#' 
#' @return
#' Return a subsetted object of the same class as \code{x}.
#' 
#' @param x object to be subsetted.
#' @param subset logical expression indicating elements or rows to keep: missing values are taken as false.
#' @param drop passed on to [ indexing operator.
#' @param reindex.all a single element logical vector indicating if residues and elements IDs have to be reindexed after subsetting.
#' @param \dots further arguments to be passed to or from other methods.
#' 
#' @seealso \code{\link[base]{subset}}, \code{\link{pdb}}, \code{\link{atoms}}, \code{\link{reindex}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' y <- subset(x, x$atoms$eleid %in% sample(x$atoms$eleid, 10))
#' is(y)
#' y <- subset(x$atoms, x$atoms$eleid %in% sample(x$atoms$eleid, 10))
#' is(y)
#' x <- coords(x)
#' y <- subset(x, x < 0)
#' is(y)
#' 
#' @keywords manip
#' 
#' @name subset.atoms
#' @export
subset.atoms <- function(x, subset, drop = FALSE, reindex.all = TRUE, ...)
{
  if (missing(subset)) 
    r <- TRUE
  else
  {
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must evaluate to logical")
    r <- r & !is.na(r)
  }
  x <- x[r, , drop = drop]
  
  if(nrow(x) == 0) x <- NULL
  
  if(reindex.all) x <- reindex.atoms(x)
  return(x)
}

#' @rdname subset.atoms
#' @export
subset.pdb <- function(x, subset, drop = FALSE, reindex.all = TRUE, ...)
{
  if (missing(subset)) 
    r <- TRUE
  else
  {
    e <- substitute(subset)
    r <- eval(e, x$atoms, parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must evaluate to logical")
    r <- r & !is.na(r)
  }
  if(any(r)) x$atoms <- x$atoms[r, , drop = drop]
  else x["atoms"] <- list(NULL)
  r <-     x$conect$eleid.1 %in% x$atoms$eleid
  r <- r & x$conect$eleid.2 %in% x$atoms$eleid
  if(any(r)) x$conect <- x$conect[r,]
  else x["conect"] <- list(NULL)
  
  if(reindex.all) x <- reindex.pdb(x)
  
  return(x)
}
