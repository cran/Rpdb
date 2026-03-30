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
#' @param crystal an object of class \sQuote{crystal} containing periodical boundary conditions used for replicating.
#' @param a.ind a vector of integers indicating the positions of the replicated cells along the a-axis.
#' @param b.ind a vector of integers indicating the positions of the replicated cells along the b-axis.
#' @param c.ind a vector of integers indicating the positions of the replicated cells along the c-axis.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}, \code{\link{crystal}}
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
replicate.coords <- function(x, crystal = NULL, a.ind = 0, b.ind = 0, c.ind = 0, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
	b = basis(x);
	if(b == "xyz") {
		check.crystal(crystal);
		x = xyz2abc(x, crystal);
	}
  
	abc.ind <- expand.grid(a.ind, b.ind, c.ind);
	x = add.abc(x, abc = abc.ind);
	if(b == "xyz") x = abc2xyz(x, crystal);
	
	return(x);
}

#' @rdname replicate
#' @export
replicate.atoms <- function(x, crystal = NULL, a.ind = 0, b.ind = 0, c.ind = 0, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
	basis = basis(x);
	if(basis == "xyz") {
		check.crystal(crystal);
		x = xyz2abc(x, crystal);
	}
	
	abc.ind = expand.grid(a.ind, b.ind, c.ind);
	
	nID = max(x$eleid);
	# 1:nrow(abc.ind)-1
	idABC = rep(seq(0, nrow(abc.ind)-1), each=natom(x));
	eleid = rep(x$eleid, nrow(abc.ind)) + idABC * nID;
	resid = rep(x$resid, nrow(abc.ind)) + idABC * max(x$resid);
  
	x = add.abc(x, abc = abc.ind);
	x$resid = resid;
	x$eleid = eleid;
	
	if(basis == "xyz") x = abc2xyz(x, crystal);
	
	return(x);
}

#' @rdname replicate
#' @export
replicate.pdb <- function(x, a.ind = 0, b.ind = 0, c.ind = 0, crystal = NULL, ...)
{
	if(! is.pdb(x)) stop("'x' must be an object of class 'pdb'");
		
	if(is.null(crystal))
		crystal <- x$crystal;
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
  na <- length(a.ind)
  nb <- length(b.ind)
  nc <- length(c.ind)
  
  ncell <- na*nb*nc
  
	eleid.1 = x$connect$eleid.1;
	eleid.2 = x$connect$eleid.2;
	nID = max(x$atoms$eleid);
	eleid.1 = rep(eleid.1, ncell) + rep(1:ncell-1, each=length(eleid.1))*nID;
	eleid.2 = rep(eleid.2, ncell) + rep(1:ncell-1, each=length(eleid.2))*nID;
	connect = connect.default(eleid.1, eleid.2);
	
	atoms = replicate.atoms(x$atoms, crystal, a.ind, b.ind, c.ind);
	
	# TODO: 1 + c();
	idABC = c(diff(range(a.ind))+1, diff(range(b.ind))+1, diff(range(c.ind))+1);
	cryst = x$crystal;
	abc   = cryst$abc * idABC;
	cryst = crystal(abc = abc, abg = cryst$abg, sgroup = cryst$sgroup);
	
	x = pdb(atoms, cryst, connect);
  
  return(x)
}

check.crystal = function(x) {
	if(is.null(x))      stop("Please specify a 'crystal' object");
	if(! is.crystal(x)) stop("'crystal' must be an object of class 'crystal'");
}

add.abc = function(x, abc) {
	L = apply(abc, 1, function(abc, x) {
			x$x1 = x$x1 + abc[1];
			x$x2 = x$x2 + abc[2];
			x$x3 = x$x3 + abc[3];
			return(x);
		}, x);
	x = do.call(rbind, L);
	return(x);
}
