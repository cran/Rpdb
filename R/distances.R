#' Inter-Atomic Distances
#' 
#' Computes inter-atomic distance vectors.
#' 
#' The purpose of the \sQuote{distances} class is to store the inter-atomic 
#' distance vectors and facilitate their manipulation when passing from the 
#' Cartesian to fractional references and vice versa.\cr The default method of 
#' the \code{distances} function is actually a builder allowing to create a 
#' \sQuote{distances} object from its different components, i.e.: \code{dx1}, 
#' \code{dx2}, \code{dx3}, and \code{basis}. All the arguments have to be 
#' specified except 'basis' which by default is set to "xyz" (Cartesian 
#' reference).
#' 
#' For objects of class \sQuote{coords}, \sQuote{atoms}, \sQuote{pdb},
#' two sets of atomic coordinates, defined by \code{sel1} and \code{sel2},
#' are extracted and inter-atomic distance vectors are computed 
#' between these two sets.
#' 
#' The method for the \sQuote{dist.point} function computes the inter-atomic distances
#' between the atoms and a specified point.
#' 
#' The method of the \code{norm} function for 
#' objects of class \sQuote{distances} computes the norm of the distances 
#' vectors. \code{type} specify how to project the distance vectors before 
#' computing the norms. By default no projection is perform. The three dx, dy, 
#' and dz components of the distance vectors are used to compute the norm. 
#' \code{type} can take the following values: \itemize{ \item   x: The distance 
#' vectors are projected over x. \item   y: The distance vectors are projected 
#' over y. \item   z: The distance vectors are projected over z. \item  xy: The 
#' distance vectors are projected in the xy-plan. \item  yz: The distance 
#' vectors are projected in the yz-plan. \item  zx: The distance vectors are 
#' projected in the zx-plan. \item xyz: The distance vectors are not projected 
#' (The three components of the distance vectors are used to compute the norm). 
#' }
#' \code{is.distances} tests if x is an object of class \sQuote{distances};
#' the test is limited to the class attribute.
#' 
#' @return The \code{distances} and \code{dist.point} functions return an object
#' of class \sQuote{distances} containing inter-atomic distance vectors.
#' The \code{norm} function return an array, with the same dimensions as the \code{dx1},
#' \code{dx2}, \code{dx3} components of the \sQuote{distances} object for which
#' the norms have to be computed, containing the norm of the distance vectors. 
#' \cr\cr
#' \code{is.distances} returns TRUE if x is an object of class \sQuote{distances}
#' and FALSE otherwise.
#' 
#' @param dx1,dx2,dx3 numeric arrays containing the first, second and third components of the distance vectors.
#' @param basis a single element character vector indicating the type of basis vector used to express the coordinates.
#' @param x,y,z an R object containing atomic coordinates.
#' @param sel1,sel2 integer or logical vectors defining two atomic selections between which the distance vectors are computed.
#' @param type a single element character vector indicating how to project the distances vectors before computing the norms. See details.
#' @param data an object of type pdb, atoms or coords, containing atomic coordinates.
#' @param subset enables sub-setting the coords object;
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso 
#' \code{\link{coords}}, \code{\link{basis}}, \code{\link{xyz2abc}}, \code{\link{abc2xyz}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' is.DCB7 <- x$atoms$resname == "DCB" & x$atoms$resid == 7
#' is.DCB8 <- x$atoms$resname == "DCB" & x$atoms$resid == 8
#' d <- distances(x, is.DCB7, is.DCB8)
#' norm(d, type = "xyz")
#' norm(d, type = "xy")
#' norm(d, type = "x")
#
#' d <- dist.point(x, c(0,0,0), subset = is.DCB7)
#' norm(d, type = "xyz")
#' 
#' @keywords classes manip
#' 
#' @name distances
#' @export
distances <- function(...)
  UseMethod("distances")

#' @rdname distances
#' @export
dist.point <- function(...)
  UseMethod("dist.point")

#' @rdname distances
#' @export
distances.default <- function(dx1 = numeric(0), dx2 = numeric(0), dx3 = numeric(0), basis = "xyz", ...) {
  if( ! is.numeric(dx1) || ! is.numeric(dx2) || ! is.numeric(dx3))
    stop("'dx1', 'dx2' and 'dx3' must be numeric");
  if(is.null(dim(dx1))) {
    if( ! is.null(dim(dx2)) && ! is.null(dim(dx3)))
      stop("'dx1', 'dx2' and 'dx3' must have the same length/dimensions");
    if(length(dx1) != length(dx2) || length(dx1) != length(dx3))
      stop("'dx1', 'dx2' and 'dx3' must have the same length");
  }
  else {
    if(any(dim(dx1) != dim(dx2)) || any(dim(dx1) != dim(dx3)))
      stop("'dx1', 'dx2' and 'dx3' must have the same dimensions");
  }
  check.basis(basis);
  
  to.return <- list(dx1 = dx1, dx2 = dx2, dx3 = dx3);
  class(to.return) <- c("distances", "list");
  attr(to.return, "basis") <- basis;
  return(to.return)
}

#' @rdname distances
#' @export
distances.coords <- function(x, sel1, sel2, ...){
  if(missing(sel1) | missing(sel2))
    stop("Please specify 'sel1' and 'sel2'")
  if(is.numeric(sel1)) {
    check.idOrLogical(sel1, natom(x), "'sel1'");
  }
  if(is.numeric(sel2)) {
    check.idOrLogical(sel2, natom(x), "'sel2'");
  }
  if(is.logical(sel1) & length(sel1) != natom(x))
    stop("'sel1' length must be equal to natom(x)")
  if(is.logical(sel2) & length(sel2) != natom(x))
    stop("'sel2' length must be equal to natom(x)")  

  xyz1 <- x[sel1,]
  xyz2 <- x[sel2,]
  
  dx1 <- t(outer(xyz2$x1, xyz1$x1, "-"))
  dx2 <- t(outer(xyz2$x2, xyz1$x2, "-"))
  dx3 <- t(outer(xyz2$x3, xyz1$x3, "-"))
  
  dimnames(dx1) <- list(sel1 = NULL, sel2 = NULL)
  dimnames(dx2) <- list(sel1 = NULL, sel2 = NULL)
  dimnames(dx3) <- list(sel1 = NULL, sel2 = NULL)
  
  to.return <- distances.default(dx1, dx2, dx3, basis = basis(x))
  return(to.return)
}

#' @rdname distances
#' @export
distances.atoms <- function(x, sel1, sel2, ...)
  distances.coords(coords(x), sel1 = sel1, sel2 = sel2, ...)

#' @rdname distances
#' @export
distances.pdb <- function(x, sel1, sel2, ...)
  distances.atoms(x$atoms, sel1 = sel1, sel2 = sel2, ...)

### Distance to a Point
#' @rdname distances
#' @export
dist.point.default <- function(data, x, y = NULL, z = NULL, subset = NULL, ...) {
	if(is.null(y)) {
		if(length(x) != 3) stop("Missing coordinates for given point!");
		y = x[2]; z = x[3]; x = x[1];
	}
	xyz = coords(data);
	if(is.null(xyz)) {
		stop("Data must be an object that contains coordinates, e.g.",
			"pdb, atoms or coords!");
	}
	basis = basis(xyz);
	if( ! is.null(subset)) {
		xyz = subset(xyz, subset);
	}
	
	dx1 = xyz$x1 - x;
	dx2 = xyz$x2 - y;
	dx3 = xyz$x3 - z;
	
	to.return <- distances.default(dx1, dx2, dx3, basis=basis);
	return(to.return)
}

#' @rdname distances
#' @export
is.distances <- function(x)
	inherits(x, "distances");

#' @rdname distances
#' @export
norm <- function(...)
  UseMethod("norm")

#' @rdname distances
#' @export
norm.distances <- function(x, type = "xyz", ...){
  if(basis(x) == "abc")
    stop("Please provide Cartesian coordinates. See 'abc2xyz'.")
  
  to.return <- switch(type,
    x   = x$dx1,
    y   = x$dx2,
    z   = x$dx3,
    xy  = sqrt(x$dx1^2 + x$dx2^2),
    yz  = sqrt(x$dx2^2 + x$dx3^2),
    zx  = sqrt(x$dx3^2 + x$dx1^2),
    xyz = sqrt(x$dx1^2 + x$dx2^2 + x$dx3^2))

  return(to.return)
}

