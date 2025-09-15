#' Projections of Objects
#' 
#' Computes 3D Projections
#' 
#' The purpose of the \sQuote{proj} helper functions is to compute
#' the spatial projections of various objects on other geometrical structures.
#' Currently, only projection of points on a line segment are implemented.
#' 
#' The methods for the \sQuote{proj.line3d} function compute the projections
#' of 1 point or a collection of points on a line segment defined by the 2
#' delimiting points.
#' 
#' The \sQuote{numeric} method projects a point, while the \sQuote{matrix} method
#' projects a collection of points on another geometric object.
#' 
#' @return The \code{proj.line3d} function returns an object
#' of class \sQuote{proj} containing the 3-dimensional coordinates of
#' the projected point, as well as a scalar value representing the fraction
#' of the path on the line segment.
#' 
#' @param p numeric array or matrix with the x,y,z coordinates of the point or collection of points.
#' @param x,y,z an R object containing the coordinates of the end-points of the line segment; can be specified either as individual values or as a matrix.
#' @param \dots currently not used.
#' 
#' @examples
#' p = c(1,2,3)
#' line = matrix(c(0,5,2,3,1,4), nrow = 2)
#' proj.line3d(p, line)
#' 
#' @keywords classes manip
#' 
#' @name proj.line3d
#' @export proj.line3d
proj.line3d <- function(p, x, y, z, ...)
  UseMethod("proj.line3d")

#' @rdname proj.line3d
#' @exportS3Method  proj.line3d numeric
proj.line3d.numeric = function(p, x, y = NULL, z = NULL, ...) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	dx = x[2] - x[1]; dx0 = p[1] - x[1];
	dy = y[2] - y[1]; dy0 = p[2] - y[1];
	dz = z[2] - z[1]; dz0 = p[3] - z[1];
	tt = (dx*dx0 + dy*dy0 + dz*dz0);
	tt = tt / (dx^2 + dy^2 + dz^2);
	t1 = 1 - tt;
	px = t1*x[1] + tt*x[2];
	py = t1*y[1] + tt*y[2];
	pz = t1*z[1] + tt*z[2];
	lst = list(x = px, y = py, z = pz, t = tt);
	class(lst) = c("proj", "list");
	return(lst);
}

#' @rdname proj.line3d
#' @exportS3Method  proj.line3d matrix
proj.line3d.matrix = function(p, x, y = NULL, z = NULL, ...) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	dx = x[2] - x[1]; dx0 = p[, 1] - x[1];
	dy = y[2] - y[1]; dy0 = p[, 2] - y[1];
	dz = z[2] - z[1]; dz0 = p[, 3] - z[1];
	tt = (dx*dx0 + dy*dy0 + dz*dz0);
	tt = tt / (dx^2 + dy^2 + dz^2);
	t1 = 1 - tt;
	px = t1*x[1] + tt*x[2];
	py = t1*y[1] + tt*y[2];
	pz = t1*z[1] + tt*z[2];
	lst = list(x = px, y = py, z = pz, t = tt);
	class(lst) = c("proj", "list");
	return(lst);
}

