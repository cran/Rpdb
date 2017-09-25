#' Helper Functions for Translation of Atomic Coordinates
#' 
#' Translation of atomic coordinates along a specific Cartesian or lattice
#' vector.
#' 
#' These functions are helper functions to perform a translation along a
#' specific Cartesian or lattice vector. All of them call either the \code{Txyz}
#' or \code{Tabc} function.
#' 
#' @return 
#' An object of the same class as \code{x} with translated coordinates.
#' 
#' @param obj an R object containing atomic coordinates.
#' @param x the x-component of the translation vector.
#' @param y the y-component of the translation vector.
#' @param z the z-component of the translation vector.
#' @param a the a-component of the translation vector.
#' @param b the b-component of the translation vector.
#' @param c the c-component of the translation vector.
#' @param mask a logical vector indicating the set of coordinates to which to apply the translation.
#' @param thickness a numeric value indicating the fraction of the thicknees of the selected atom to be added to the translation vector (Usually 0, 0.5 or 1. See details).
#' @param cryst1 an object of class \sQuote{cryst1} use to convert Cartesian into fraction coordinates (or Vis Versa) when need.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso 
#' \code{\link{Txyz}}, \code{\link{Tabc}}\cr
#' Passing from Cartesian to fractional coordinates (or Vis Versa):\cr
#' \code{\link{xyz2abc}}, \code{\link{abc2xyz}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#' visualize(x, mode = NULL)
#' visualize(Ty(x, 10), mode = NULL)
#' visualize(Ty(x, 10, mask=x$atoms$resid==1), mode = NULL)
#' visualize(Tb(x, 1 ), mode = NULL)
#' visualize(Tb(x, 1 , mask=x$atoms$resid==1), mode = NULL)
#' 
#' # Lets build a C70/Pentacene dimer with an inter-molecular distance equal to 3.5
#' C70 <- read.pdb(system.file("examples/C70.pdb",package="Rpdb"))
#' Pen <- read.pdb(system.file("examples/Pentacene.pdb",package="Rpdb"))
#' x <- merge(C70,Pen)
#' visualize(x, mode = NULL)
#' viewXY()
#' visualize(Tz(x, z=3.5, mask=x$atoms$resname=="C70", thickness=0.5), mode = NULL)
#' viewXY()
#' 
#' @keywords manip
#' 
#' @name translationHelpers
#' @export
Tx <- function(...)
  UseMethod("Tx")
#' @rdname translationHelpers
#' @export
Tx.coords <- function(obj, x = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...)
  Txyz(obj, x = x, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
#' @rdname translationHelpers
#' @export
Tx.pdb <- function(obj, x = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...)
  Txyz(obj, x = x, mask = mask, thickness = thickness, cryst1 = cryst1, ...)

#' @rdname translationHelpers
#' @export
Ty <- function(...)
  UseMethod("Ty")
#' @rdname translationHelpers
#' @export
Ty.coords <- function(obj, y = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...)
  Txyz(obj, y = y, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
#' @rdname translationHelpers
#' @export
Ty.pdb <- function(obj, y = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...)
  Txyz(obj, y = y, mask = mask, thickness = thickness, cryst1 = cryst1, ...)

#' @rdname translationHelpers
#' @export
Tz <- function(...)
  UseMethod("Tz")
#' @rdname translationHelpers
#' @export
Tz.coords <- function(obj, z = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...)
  Txyz.coords(obj, z = z, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
#' @rdname translationHelpers
#' @export
Tz.pdb <- function(obj, z = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...)
  Txyz.pdb(obj, z = z, mask = mask, thickness = thickness, cryst1 = cryst1, ...)

#' @rdname translationHelpers
#' @export
Ta <- function(...)
  UseMethod("Ta")
#' @rdname translationHelpers
#' @export
Ta.coords <- function(obj, a = 0, mask = TRUE, cryst1 = NULL, ...)
  Tabc(obj, a = a, mask = mask, cryst1 = cryst1, ...)
#' @rdname translationHelpers
#' @export
Ta.pdb <- function(obj, a = 0, mask = TRUE, cryst1 = obj$cryst1, ...)
  Tabc(obj, a = a, mask = mask, cryst1 = cryst1, ...)

#' @rdname translationHelpers
#' @export
Tb <- function(...)
  UseMethod("Tb")
#' @rdname translationHelpers
#' @export
Tb.coords <- function(obj, b = 0, mask = TRUE, cryst1 = NULL, ...)
  Tabc(obj, b = b, mask = mask, cryst1 = cryst1, ...)
#' @rdname translationHelpers
#' @export
Tb.pdb <- function(obj, b = 0, mask = TRUE, cryst1 = obj$cryst1, ...)
  Tabc(obj, b = b, mask = mask, cryst1 = cryst1, ...)

#' @rdname translationHelpers
#' @export
Tc <- function(...)
  UseMethod("Tc")
#' @rdname translationHelpers
#' @export
Tc.coords <- function(obj, c = 0, mask = TRUE, cryst1 = NULL, ...)
  Tabc.coords(obj, c = c, mask = mask, cryst1 = cryst1, ...)
#' @rdname translationHelpers
#' @export
Tc.pdb <- function(obj, c = 0, mask = TRUE, cryst1 = obj$cryst1, ...)
  Tabc.pdb(obj, c = c, mask = mask, cryst1 = cryst1, ...)

