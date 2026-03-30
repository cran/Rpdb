#' Graphics Tools
#' 
#' Tools to modify or set \code{rgl} objects, scenes and the observer.


#' \code{translate} performs a translation of the observer to a new position
#'   specified using relative shifts.

#' \code{trans} is an abbreviated alias to function \code{translate}.

#' @return \code{translate} returns invisibly the old position of the observer.

#' @param x,y,z amount to shift the observer in the x, y and z directions.
#    Negative shifts are to the left, bottom and out of the scene,
#    while positive shifts are to the right, top and into the scene.
#' @param auto If TRUE, the location will be set automatically by RGL
#'   to make the whole bounding box visible.
#'   Unfortunately, it may fail, as it seems there are some bugs in rgl.


### Translate

#' @name translate
#' @export
trans = function(x, y = 0, z = 0, auto = FALSE) {
	old = translate(x, y, z, auto=auto);
}

#' @name translate
#' @export
translate = function(x, y = 0, z = 0, auto = FALSE) {
	old = observer3d();
	observer3d(old + c(x, y, z), auto=auto);
	invisible(old);
}
