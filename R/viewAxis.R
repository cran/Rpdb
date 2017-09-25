#' Set the View of the \sQuote{rgl} Scene
#' 
#' Set the view of the current \sQuote{rgl} scene aligning one vector 
#' perpendicularly to the screen and placing another in the horizontal plan.
#' 
#' \code{viewAxis} set the view of the current rgl scene (by setting
#' \code{UserMatrix}. See \code{\link[rgl]{par3d}} for more details) so that
#' \code{V1} is perpendicular to the screen and \code{V2} is in the horizontal
#' plan. The other functions documented here are helper functions calling
#' \code{viewAxis} to set the view using particular Cartesian or lattice
#' vectors. For functions \code{viewAB}, \code{viewBC} and \code{viewCA} a
#' \sQuote{cryst1} object has to be specifyed to defined the lattice vectors
#' used to set the view. The function \code{viewInertia} computes the inertia
#' tensor from atomic coordinates and masses (see \code{\link{inertia}}) and set
#' the view to its eigen vectors basis set.
#' 
#' @param V1 a length 3 numeric vector.
#' @param V2 a length 3 numeric vector.
#' @param cryst1 an object of class \sQuote{cryst1}.
#' @param x an R object containing atomic coordinates.
#' @param m a numeric vector containing atomic masses.
#' 
#' @seealso \code{\link{visualize}}, \code{\link{cell.coords}}, \code{\link[rgl]{par3d}}, \code{\link[rgl]{rgl.open}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"))
#' visualize(x, mode = NULL)
#' viewAB(x$cryst1)
#' 
#' C70 <- read.pdb(system.file("examples/C70.pdb",package="Rpdb"))
#' visualize(C70, mode = NULL)
#' viewXY()
#' viewInertia(C70)
#' 
#' @keywords dynamic
#' 
#' @name viewAxis
#' @export
viewAxis <- function(V1, V2){
  if(is.list(V1)|is.list(V2))
    stop("'V1' and 'V2' can not be lists")
  if(length(V1)!=3|length(V2)!=3)
    stop("'V1' and 'V2' must be of length 3")
  if(vectNorm(V1)==0|vectNorm(V2)==0)
    stop("'V1' and 'V2' can not be null vectors")
  V1 <- matrix(c(V1,0))
  V2 <- matrix(c(V2,0))
  A <- acos(V1[1]/sqrt(sum(V1[c(1,2)]^2)))
  A <- ifelse(V1[2] > 0, -A, A)
  R1 <- rgl::rotationMatrix(A ,0 , 0, 1)
  V1 <- R1%*%V1
  V2 <- R1%*%V2
  A <- acos(V1[3]/sqrt(sum(V1[c(1,3)]^2)))
  A <- ifelse(V1[1] > 0, -A, A)
  R2 <- rotationMatrix(A ,0 , 1, 0)
  Vz <- R2%*%V1
  V2 <- R2%*%V2
  A <- acos(V2[1]/sqrt(sum(V2[c(1,2)]^2)))
  A <- ifelse(V2[2] > 0, -A, A)
  R3 <- rgl::rotationMatrix(A ,0 , 0, 1)
  rgl::par3d(userMatrix=R3%*%R2%*%R1)
}

#' @rdname viewAxis
#' @export
viewXY <- function()
  viewAxis(c(1,0,0),c(0,1,0))

#' @rdname viewAxis
#' @export
viewYZ <- function()
  viewAxis(c(0,1,0),c(0,0,1))

#' @rdname viewAxis
#' @export
viewZX <- function()
  viewAxis(c(0,0,1),c(1,0,0))

#' @rdname viewAxis
#' @export
viewAB <- function(cryst1){
  if(missing(cryst1))
    stop("Please specify 'cryst1' to defined the lattice vectors used to set the view")
  cell <- cell.coords.cryst1(cryst1)
  viewAxis(cell[,"a"],cell[,"b"])
}

#' @rdname viewAxis
#' @export
viewBC <- function(cryst1){
  if(missing(cryst1))
    stop("Please specify 'cryst1' to defined the lattice vectors used to set the view")
  cell <- cell.coords.cryst1(cryst1)
  viewAxis(cell[,"b"],cell[,"c"])
}

#' @rdname viewAxis
#' @export
viewCA <- function(cryst1){
  if(missing(cryst1))
    stop("Please specify 'cryst1' to defined the lattice vectors used to set the view")
  cell <- cell.coords.cryst1(cryst1)
  viewAxis(cell[,"c"],cell[,"a"])
}

#' @rdname viewAxis
#' @export
viewInertia <- function(x, m = NULL){
  M <- diag(4)
  M[1:3,1:3] <- eigen(inertia(x, m))$vectors
  rgl::par3d(userMatrix=M)
}
