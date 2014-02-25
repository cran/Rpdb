#  Merge by row two data frames of class 'coords'.

merge.coords <- function(x, y, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.coords(y)) stop("'y' must be an object of class 'coords'")
  
  if(basis(x) != basis(y)) stop("'x' and 'y' basis differ")
  
  to.return <- rbind(x,y)
  return(to.return)
}

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
