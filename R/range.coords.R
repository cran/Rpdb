#' Range of Atomic Coordinates
#' 
#' Determines the range of atomic coordinates.
#' 
#' \code{range} computes the range of each of the coordinates.
#' \code{lattice.range} computes the range of the relative coordinates when expressed
#'   in units of lattice vectors.
#' 
#' @return Return a \code{\link{data.frame}} whose columns contain the range of
#'   the first, second and third coordinates of \code{x}.
#'   
#' @param x an R object containing atomic coordinates.
#' @param xyz set of lattice vectors, by default extracted from the PDB object \code{x}.
#' @param na.rm logical, indicating if \code{\link{NA}}'s should be omitted.
#' @param finite logical, indicating if all non-finite elements should be omitted.
#' @param \dots further arguments passed to or from other methods. 
#'
#' @seealso \code{range}, \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' range(x)
#' range(range(x))
#' 
#' @keywords manip


### Functions:

#' @name range.coords
#' @export
range.coords <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
	if(! is.coords(x)) stop("'x' must be an object of class 'coords'");
	
	to.return = lapply(x[, c("x1","x2","x3")], range, na.rm, finite);
	to.return = as.data.frame(to.return, row.names = c("min","max"));

  c.names <- unlist(strsplit(basis(x), ""))
  colnames(to.return) <- c.names
  
  return(to.return)
}

#' @rdname range.coords
#' @export
range.atoms <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  to.return <- range.coords(coords(x), na.rm, finite)
  return(to.return)
}

#' @rdname range.coords
#' @export
range.pdb = function(x, na.rm = FALSE, finite = FALSE, ...)
{
	if(! is.pdb(x)) stop("'x' must be an object of class 'pdb'");
	
	res = range.coords(coords(x), na.rm=na.rm, finite=finite);
	return(res);
}


### Range of Coordinates as Lattice Vectors

# Hack to enable methods(lattice);
lattice = function() {}

#' @name range.coords
#' @export lattice.range
lattice.range = function(x, xyz = NULL, ..., na.rm = TRUE)
	UseMethod("lattice.range");

#' @rdname range.coords
#' @exportS3Method  lattice.range pdb
lattice.range.pdb = function(x, xyz = NULL, ..., na.rm = TRUE) {
	range.lattice.pdb(x, xyz=xyz, ..., na.rm=na.rm);
}


#' @rdname range.coords
#' @method range lattice
#' @export range.lattice
range.lattice = function(x, ...) {
	UseMethod("range.lattice");
}

#' @rdname range.coords
#' @method range.lattice pdb
#' @exportS3Method range.lattice pdb
range.lattice.pdb = function(x, xyz = NULL, ..., na.rm = TRUE) {
	if(is.null(xyz)) {
		cell = get.PDBCrystal(x);
		xyzV = cell.coords(cell)
	} else xyzV = xyz;
	xyzE = coords(x$atoms);
	tt = solve(xyzV, t(xyzE));
	tt = apply(tt, 1, range, na.rm=na.rm);
	# xyz = xyzV %*% t(tt);
	return(tt);
}
