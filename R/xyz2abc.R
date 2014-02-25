#  Convert Cartesian coordinates into fractional atoms and vis versa.

xyz2abc <- function(...)
  UseMethod("xyz2abc")

xyz2abc.coords <- function(x, cryst1, ...)
{
  if(missing(cryst1)) stop("Please specify a 'cryst1' object")
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  if(basis(x) != "xyz") stop("Coordinates are not 'xyz' atoms")
  
  cell <- cell.coords(cryst1)
  to.return <- solve(cell)%*%t(x)
  to.return <- coords.default(to.return["a",],to.return["b",],to.return["c",],"abc")

  return(to.return)
  
}

xyz2abc.atoms <- function(x, cryst1, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  value <- xyz2abc.coords(coords.atoms(x), cryst1)
  coords(x) <- value
  
  return(x)
}

xyz2abc.pdb <- function(x, cryst1 = x$cryst1, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  value <- xyz2abc.coords(coords.pdb(x), cryst1)
  coords(x) <- value
  
  return(x)
}

abc2xyz <- function(...)
  UseMethod("abc2xyz")

abc2xyz.coords <- function(x, cryst1, ...)
{
  if(missing(cryst1)) stop("Please specify a 'cryst1' object")
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  if(basis.default(x) != "abc") stop("Coordinates are not 'abc' atoms")
  
  cell <- cell.coords(cryst1)
  to.return <- cell%*%t(x)
  to.return <- coords.default(to.return["x",],to.return["y",],to.return["z",],"xyz")
  
  return(to.return)
  
}

abc2xyz.atoms <- function(x, cryst1, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  value <- abc2xyz.coords(coords.atoms(x), cryst1)
  coords(x) <- value
  
  return(x)
}

abc2xyz.pdb <- function(x, cryst1 = x$cryst1, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  value <- abc2xyz.coords(coords.pdb(x), cryst1)
  coords(x) <- value
  return(x)
}

