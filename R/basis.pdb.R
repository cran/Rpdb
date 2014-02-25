#  Functions to get or set the basis of an object of class 'pdb'.

basis.pdb <- function(x)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  to.return <- basis.default(x$atoms)
  return(to.return)
}

'basis<-.pdb' <- function(x, value)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  if(!value %in% c("xyz","abc")) stop("Unrecognized basis set")
  attr(x$atoms, "basis") <- value
  return(x)
}
