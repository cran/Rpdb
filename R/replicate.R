#' Replicate Atomic Coordinates
#' 
#' Replicate atomic coordinates using periodic boundary conditions.
#' 
#' The \code{replicate} function replicate a unit cell along the lattice vectors a, b and c
#' as as many times as indicated by the \code{a.ind}, \code{b.ind} and \code{c.ind} arguments.
#' Discontinuous integer vectors can be used for \code{a.ind}, \code{b.ind} and \code{c.ind}
#' to create layered supercells (See examples).
#' 
#' @return Return an object of class \sQuote{pdb} with replicated atomic coordinates.
#' 
#' @param x an R object containing atomic coordinates to be replicated.
#' @param cryst1 an object of class \sQuote{cryst1} containing periodical boundary conditions used for replicating.
#' @param a.ind a vector of integers indicating the positions of the replicated cells along the a-axis.
#' @param b.ind a vector of integers indicating the positions of the replicated cells along the b-axis.
#' @param c.ind a vector of integers indicating the positions of the replicated cells along the c-axis.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}, \code{\link{cryst1}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' 
#' # Create a 3x3 supercell
#' y <- replicate(x, a.ind = 0:2, b.ind = 0:2, c.ind = 0:2)
#' 
#' # Create a 3x3 supercell which might need to be wrapped (if some molecules are outside the cell)
#' y <- replicate(x, a.ind = -1:1, b.ind = -1:1, c.ind = -1:1)
#' 
#' # Create a layered supercell with a vacuum layer in the bc-plan
#' y <- replicate(x, a.ind = c(0,2), b.ind = 0:2, c.ind = 0:2)
#'
#' @keywords manip
#' 
#' @name replicate
#' @export
replicate <- function(x, ...)
  UseMethod("replicate")

#' @rdname replicate
#' @export
replicate.coords <- function(x, cryst1 = NULL, a.ind = 0, b.ind = 0, c.ind = 0, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
  b <- basis(x)
  if(b == "xyz")
  {
    if( is.null(cryst1))   stop("Please specify a 'cryst1' object")
    if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
    x <- xyz2abc(x, cryst1) 
  }
  
  abc.ind <- expand.grid(a.ind, b.ind, c.ind)
  
  L <- apply(abc.ind, 1,
             function(abc, x)
             {
               x$x1 <- x$x1 + abc[1]
               x$x2 <- x$x2 + abc[2]
               x$x3 <- x$x3 + abc[3]
               return(x)
             }, x)
  
  x <- do.call(rbind, L)
  if(b == "xyz") x <- abc2xyz(x, cryst1)
  
  return(x)
}

#' @rdname replicate
#' @export
replicate.atoms <- function(x, cryst1 = NULL, a.ind = 0, b.ind = 0, c.ind = 0, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
  basis <- basis(x)
  if(basis == "xyz")
  {
    if( is.null(cryst1))   stop("Please specify a 'cryst1' object")
    if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
    x <- xyz2abc(x, cryst1) 
  }
  
  abc.ind <- expand.grid(a.ind, b.ind, c.ind)
  
  L <- apply(abc.ind, 1,
             function(abc, x)
             {
               x$x1 <- x$x1 + abc[1]
               x$x2 <- x$x2 + abc[2]
               x$x3 <- x$x3 + abc[3]
               return(x)
             }, x)
  
  eleid <- rep(x$eleid, nrow(abc.ind)) + rep(1:nrow(abc.ind)-1,each=natom(x))*max(x$eleid)
  resid <- rep(x$resid, nrow(abc.ind)) + rep(1:nrow(abc.ind)-1,each=natom(x))*max(x$resid)
  
  x <- do.call(rbind, L)
  x$resid <- resid
  x$eleid <- eleid
  
  if(basis == "xyz") x <- abc2xyz(x, cryst1)
  
  return(x)
}

#' @rdname replicate
#' @export
replicate.pdb <- function(x, a.ind = 0, b.ind = 0, c.ind = 0, cryst1 = NULL, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  if(is.null(cryst1))
    cryst1 <- x$cryst1
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
  na <- length(a.ind)
  nb <- length(b.ind)
  nc <- length(c.ind)
  
  ncell <- na*nb*nc
  
  eleid.1 <- x$conect$eleid.1
  eleid.2 <- x$conect$eleid.2
  eleid.1 <- rep(eleid.1, ncell) + rep(1:ncell-1,each=length(eleid.1))*max(x$atoms$eleid)
  eleid.2 <- rep(eleid.2, ncell) + rep(1:ncell-1,each=length(eleid.2))*max(x$atoms$eleid)
  conect <- conect.default(eleid.1, eleid.2)
  
  atoms <- replicate.atoms(x$atoms, cryst1, a.ind, b.ind, c.ind)
  
  abc <- x$cryst1$abc*c(diff(range(a.ind))+1, diff(range(b.ind))+1, diff(range(c.ind))+1)
  cryst1 <- cryst1(abc = abc, abg = x$cryst1$abg, sgroup = x$cryst1$sgroup)

  x <- pdb(atoms, cryst1, conect)
  
  return(x)
}