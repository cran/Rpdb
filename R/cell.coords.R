#  Compute the Cartesian coordinates of lattice vectors.

cell.coords <- function(...)
  UseMethod("cell.coords")

cell.coords.default <- function(abc, abg = c(90,90,90), digits = 3, ...)
{
  if(missing(abc)) stop("Please provide at list a 'abc' vector containing the length of the lattice vectors")
  if(length(abc) != 3) stop("'abc' must be a vector of length 3")
  if(length(abg) != 3) stop("'abg' must be a vector of length 3")
  
  abg <- abg*pi/180
  
  M <- matrix(ncol=3,nrow=3)
  M[ ,1] <- c(abc[1],0,0)
  M[ ,2] <- c(abc[2]*cos(abg[3]),abc[2]*sin(abg[3]),0)
  M[1,3] <-   abc[3]*(cos(abg[2]))
  M[2,3] <-   abc[3]*(cos(abg[1])-cos(abg[2])*cos(abg[3]))/sin(abg[3])
  M[3,3] <-   abc[3]*sqrt(1+2*cos(abg[1])*cos(abg[2])*cos(abg[3])-(cos(abg[1]))^2-(cos(abg[2]))^2-(cos(abg[3]))^2)/sin(abg[3])
  M <- round(M, digits=3)
  dimnames(M) <- list(c("x","y","z"), c("a","b","c"))
  
  return(M)
}

cell.coords.cryst1 <- function(x, digits = 3, ...)
{
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
  
  M <- cell.coords.default(x$abc, x$abg, digits)
  return(M)
}

cell.coords.pdb <- function(x, digits = 3, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'atoms'")
  if(is.null(x$cryst1)) stop("'x' must contained a 'cryst1' object")
  
  M <- cell.coords.cryst1(x$cryst1, digits = 3)
  return(M)
}

