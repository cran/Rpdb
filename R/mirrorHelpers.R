#' Helper Functions for reflection of Atomic Coordinates
#' 
#' Reflection of atomic coordinates with respect to a specific Cartesian plane or
#' a plane defined by two lattice vectors.
#' 
#' These functions are helper functions to perform a reflection with respect to 
#' a specific Cartesian plan or a plan defined by two lattice vectors. All of 
#' them call the \code{mirror} function.
#' 
#' @return An object of the same class as \code{x} with reflected coordinates.
#'   
#' @param x an R object containing atomic coordinates.
#' @param mask a logical vector indicating the set of coordinates to which to
#'   apply the reflection.
#' @param crystal an object of class \sQuote{crystal} used to convert fractional
#'   into Cartesian coordinates (when needed).
#' @param \dots further arguments passed to or from other methods.
#'   
#' @seealso \code{\link{mirror}} and \code{\link{xyz2abc}},
#' \code{\link{abc2xyz}} for passing from Cartesian to fractional coordinates
#' (or Vice Versa).
#' 
#' @examples
#' # First lets read a pdb file
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' visualize(x, mode = NULL)
#'
#' # Mirror operation with respect to the ab-plane
#' visualize(Mab(x), mode = NULL)
#' # Mirror operation with respect to the ab-plane for residue 1
#' visualize(Mab(x, mask = x$atoms$resid == 1), mode = NULL)
#' 
#' @keywords manip
#' 
#' @name mirrorHelpers
#' @export
Mxy <- function(...)
  UseMethod("Mxy")

#' @rdname mirrorHelpers
#' @export
Mxy.coords = function(x, mask = TRUE, crystal = NULL, ...)
	mirror(x, c(0,0,0), c(1,0,0), c(0,1,0), mask=mask, crystal=crystal);

#' @rdname mirrorHelpers
#' @export
Mxy.pdb = function(x, mask = TRUE, crystal = x$crystal, ...)
	mirror(x, c(0,0,0), c(1,0,0), c(0,1,0), mask=mask, crystal=crystal);


#' @rdname mirrorHelpers
#' @export
Myz <- function(...)
  UseMethod("Myz")

#' @rdname mirrorHelpers
#' @export
Myz.coords = function(x, mask = TRUE, crystal = NULL, ...)
	mirror(x, c(0,0,0), c(0,1,0), c(0,0,1), mask=mask, crystal=crystal);

#' @rdname mirrorHelpers
#' @export
Myz.pdb = function(x, mask = TRUE, crystal = x$crystal, ...)
	mirror(x, c(0,0,0), c(0,1,0), c(0,0,1), mask=mask, crystal=crystal);


#' @rdname mirrorHelpers
#' @export
Mzx <- function(...)
  UseMethod("Mzx")

#' @rdname mirrorHelpers
#' @export
Mzx.coords = function(x, mask = TRUE, crystal = NULL, ...)
	mirror(x, c(0,0,0), c(0,0,1), c(1,0,0), mask=mask, crystal=crystal);

#' @rdname mirrorHelpers
#' @export
Mzx.pdb = function(x, mask = TRUE, crystal = x$crystal, ...)
	mirror(x, c(0,0,0), c(0,0,1), c(1,0,0), mask=mask, crystal=crystal);


#' @rdname mirrorHelpers
#' @export
Mab <- function(...)
  UseMethod("Mab")

#' @rdname mirrorHelpers
#' @export
Mab.coords = function(x, crystal, mask = TRUE, ...) {
	cell = get.cellMirrorMiss(x, crystal=crystal);
	mirror(x, c(0,0,0), cell[,"a"], cell[,"b"], mask=mask, crystal=crystal);
}
  
#' @rdname mirrorHelpers
#' @export
Mab.pdb = function(x, crystal = x$crystal, mask = TRUE, ...) {
	cell = get.cellMirror(x, crystal=crystal);
	mirror(x, c(0,0,0), cell[,"a"], cell[,"b"], mask=mask, crystal=crystal);
}

#' @rdname mirrorHelpers
#' @export
Mbc <- function(...)
  UseMethod("Mbc")

#' @rdname mirrorHelpers
#' @export
Mbc.coords = function(x, crystal, mask = TRUE, ...) {
	cell = get.cellMirrorMiss(x, crystal=crystal);
	mirror(x, c(0,0,0), cell[,"b"], cell[,"c"], mask=mask, crystal=crystal);
}

#' @rdname mirrorHelpers
#' @export
Mbc.pdb <- function(x, crystal = x$crystal, mask = TRUE, ...) {
	cell = get.cellMirror(x, crystal=crystal);
	mirror(x, c(0,0,0), cell[,"b"], cell[,"c"], mask=mask, crystal=crystal);
}

#' @rdname mirrorHelpers
#' @export
Mca <- function(...)
  UseMethod("Mca")

#' @rdname mirrorHelpers
#' @export
Mca.coords = function(x, crystal, mask = TRUE, ...) {
	cell = get.cellMirrorMiss(x, crystal=crystal);
	mirror(x, c(0,0,0), cell[,"c"], cell[,"a"], mask=mask, crystal=crystal);
}

#' @rdname mirrorHelpers
#' @export
Mca.pdb = function(x, crystal = x$crystal, mask = TRUE, ...) {
	cell = get.cellMirror(x, crystal=crystal);
	mirror(x, c(0,0,0), cell[,"c"], cell[,"a"], mask=mask, crystal=crystal);
}

### Helper

# Cell:
get.cellMirror = function(x, crystal) {
	check.crystalMirror(crystal);
	cell.coords(x, crystal=crystal);
}
get.cellMirrorMiss = function(x, crystal) {
	check.crystalMirrorMiss(crystal);
	cell.coords(x, crystal=crystal);
}

# Check:
check.crystalMirror = function(x) {
	if(is.null(x))
		stop("'crystal' is required to define the mirror plane!");
}
check.crystalMirrorMiss = function(x) {
	if(missing(x) || is.null(x))
		stop("'crystal' is required to define the mirror plane!");
}
