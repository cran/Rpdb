#' Add Axes or PBC Box to the \sQuote{rgl} Scene
#' 
#' Add lattice vectors, Cartesian axes or PBC box to the current \sQuote{rgl} 
#' scene.
#' 
#' \code{addABC}: Add the lattice vectors a, b and c to the current rgl 
#' device.\cr \code{addXYZ}: Add the Cartesian axes x, y and z to the current 
#' rgl device.\cr \code{addPBCBox}: Add a box representing the Periodic Boundary
#' Conditions of a molecular system.
#' 
#' @return Return (using invisible) a two-column data.frame containing the IDs 
#'   and type indicators of the objects added to the scene.
#'   
#' @param x an object of class \sQuote{cryst1} containing unit cell parameters.
#' @param lwd a numeric value indicating the line width used to draw the axes or
#'   the PBC box.
#' @param labels a logical value indicating whether the labels of the axes have 
#'   to be drawn.
#' @param cex a numeric value indicating the magnification used to draw the 
#'   labels of the axes.
#'   
#' @seealso \code{\link{visualize}}, \code{\link[rgl]{rgl.open}}, \code{\link[rgl]{par3d}},
#'   \code{\link{addLabels}}
#'   
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' visualize(x, type = "l", xyz = FALSE, abc = FALSE, pbc.box = FALSE, mode = NULL)
#' addXYZ()
#' addABC(x$cryst1)
#' addPBCBox(x$cryst1)
#' 
#' @keywords dynamic
#'     
#' @name addAxes
#' @export
addABC <- function(x, lwd = 2, labels = TRUE, cex = 2){
  if(missing(x)) stop("Please specify a 'cryst1' object")
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1")
  
  cell <- cell.coords(x)
  
  seg.id <- rgl::segments3d(
    rbind(
      c(0,0,0),cell[,1],
      c(0,0,0),cell[,2],
      c(0,0,0),cell[,3]
    ),
    col=c("red","red","green","green","blue","blue"),
    lwd=lwd
  )
  seg.id <- data.frame(id = seg.id, type = "abc.seg")
  
  an <- cell[,1]/sqrt(sum(cell[,1]^2))
  bn <- cell[,2]/sqrt(sum(cell[,2]^2))
  cn <- cell[,3]/sqrt(sum(cell[,3]^2))
  
  lab.id <- NULL
  if(labels){
    lab.id <- rgl::text3d(
      rbind(
        cell[,1]+an*1.2,
        cell[,2]+bn*1.2,
        cell[,3]+cn*1.2
      ),
      texts = c("a","b","c"),
      col = c("red","green","blue"),
      cex=cex
    )
    lab.id <- data.frame(id = lab.id, type = "abc.lab")
  }
  ids <- rbind(seg.id, lab.id)
  
  invisible(ids)
}

#' @rdname addAxes
#' @export
addXYZ <- function(lwd = 2, labels= TRUE, cex = 2){
  seg.id <- rgl::segments3d(
    rbind(
      c(0,0,0), c( 5, 0, 0),
      c(0,0,0), c( 0, 5, 0),
      c(0,0,0), c( 0, 0, 5),
      c(0,0,0), c(-5, 0, 0),
      c(0,0,0), c( 0,-5, 0),
      c(0,0,0), c( 0, 0,-5)
    ),
    lwd = lwd,
    alpha=c(rep(1, 6), rep(0, 6))
  )
  seg.id <-data.frame(id = seg.id, type = "xyz.seg")
  
  lab.id <- NULL
  if(labels){
    lab.id <- rgl::text3d(
      rbind(
        c(0,0,0) + c(5.0+1.0,0.0    ,0.0    ),
        c(0,0,0) + c(0.0    ,5.0+1.0,0.0    ),
        c(0,0,0) + c(0.0    ,0.0    ,5.0+1.0)
      ),
      texts=c("x","y","z"),
      cex = cex
    )
    lab.id <- data.frame(id = lab.id, type = "xyz.lab")
  }
  ids <- rbind(seg.id, lab.id)
  
  invisible(ids)
}

#' @rdname addAxes
#' @export
addPBCBox <- function(x, lwd = 2){
  if(missing(x)) stop("Please specify a 'cryst1' object")
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
  
  cell <- cell.coords(x)
  
  seg.id <- rgl::segments3d(
    rbind(
      c(0,0,0)                  , cell[,1],
      c(0,0,0)                  , cell[,2],
      c(0,0,0)                  , cell[,3],
      cell[,1]+cell[,2]         , cell[,1],
      cell[,1]+cell[,2]         , cell[,2],
      cell[,1]+cell[,2]         , cell[,1]+cell[,2]+cell[,3],
      cell[,1]+cell[,3]         , cell[,1],cell[,1]+cell[,3],
      cell[,1]+cell[,3]+cell[,2], cell[,1]+cell[,3],cell[,3],
      cell[,2]+cell[,3]         , cell[,2]+cell[,3]+cell[,1],
      cell[,2]+cell[,3]         , cell[,2],
      cell[,2]+cell[,3]         , cell[,3]
    ),
    col = "black",
    lwd = lwd
  )
  seg.id <- data.frame(id = seg.id, type = "pbc.box")
  
  invisible(seg.id)
}
