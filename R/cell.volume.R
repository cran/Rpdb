cell.volume <- function(...)
  UseMethod("cell.volume")

cell.volume.cryst1 <- function(x, ...)
{
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
   
  V <-  prod(x$abc)*sqrt(1 - sum(cos(x$abg)^2) + 2*prod(cos(x$abg)))
#   attr(V, "unit") <- "AngtromCube"
  return(V)
}

cell.volume.pdb <- function(x, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'atoms'")
  if(is.null(x$cryst1)) stop("'x' must contained a 'cryst1' object")
  
  V <- cell.volume.cryst1(x$cryst1)
  return(V)
}
