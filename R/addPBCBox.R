addPBCBox <- function(x, lwd = 2){
  if(missing(x)) stop("Please specify a 'cryst1' object")
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
  
  cell <- cell.coords(x)
  
  seg.id <- segments3d(
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