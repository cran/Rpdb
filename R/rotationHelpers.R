#' Helper Functions for Rotation of Atomic Coordinates
#' 
#' Rotation of atomic coordinates along a specific Cartesian vector.
#' 
#' These functions are helper functions to perform a rotation around a specific 
#' Cartesian vector. All of them call the \code{R} function.
#' 
#' @return An object of the same class as \code{x} with rotated coordinates.
#'   
#' @param x an R object containing atomic coordinates.
#' @param angle the angle of the rotation in degrees.
#' @param mask a logical vector indicating the set of coordinates to which the
#'   rotation has to be applyed.
#' @param cryst1 an object of class \sQuote{cryst1} use to convert fractional
#'   into Cartesian coordinates when need.
#' @param \dots further arguments passed to or from other methods.
#'   
#' @seealso \code{\link{R}} and \code{\link{xyz2abc}}, \code{\link{abc2xyz}} for
#'   passing from Cartesian to fractional coordinates (or Vis Versa).
#' 
#' @examples 
#' # First lets read a pdb file
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' cell <- cell.coords(x)
#' visualize(x, mode = NULL)
#'
#' # Rotation of the structure around the z-axis
#' visualize(Rz(x, 90), mode = NULL)
#' # Rotation of the residue 1 around the c-axis
#' visualize(Rz(x, 90, mask = x$atoms$resid == 1), mode = NULL)
#' 
#' @keywords manip
#' 
#' @name rotationHelpers
#' @export
Rx <- function(...)
  UseMethod("Rx")

#' @rdname rotationHelpers
#' @export
Rx.coords <- function(x, angle = 0, mask = TRUE, cryst1 = NULL, ...)
  R(x, angle = angle, x = 1, y = 0, z = 0, mask = mask, cryst1 = cryst1, ...)

#' @rdname rotationHelpers
#' @export
Rx.pdb <- function(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, ...)
  R(x, angle = angle, x = 1, y = 0, z = 0, mask = mask, cryst1 = cryst1, ...)


#' @rdname rotationHelpers
#' @export
Ry <- function(...)
  UseMethod("Ry")

#' @rdname rotationHelpers
#' @export
Ry.coords <- function(x, angle = 0, mask = TRUE, cryst1 = NULL, ...)
  R(x, angle = angle, x = 0, y = 1, z = 0, mask = mask, cryst1 = cryst1, ...)

#' @rdname rotationHelpers
#' @export
Ry.pdb <- function(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, ...)
  R(x, angle = angle, x = 0, y = 1, z = 0, mask = mask, cryst1 = cryst1, ...)


#' @rdname rotationHelpers
#' @export
Rz <- function(...)
  UseMethod("Rz")

#' @rdname rotationHelpers
#' @export
Rz.coords <- function(x, angle = 0, mask = TRUE, cryst1 = NULL, ...)
  R(x, angle = angle, x = 0, y = 0, z = 1, mask = mask, cryst1 = cryst1, ...)

#' @rdname rotationHelpers
#' @export
Rz.pdb <- function(x, angle = 0, mask = TRUE, cryst1 = x$cryst1, ...)
  R(x, angle = angle, x = 0, y = 0, z = 1, mask = mask, cryst1 = cryst1, ...)
