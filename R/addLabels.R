#' Add Labels to the \sQuote{rgl} Scene
#' 
#' Add Labels to the current \sQuote{rgl} scene.
#' 
#' \code{addResLab} add residue labels to the scene. If \code{at.centre==TRUE} 
#' only one label per residue is added at the centre of the residue. Otherwise, 
#' residue labels are added at each atomic positions. \code{addEleLab} add 
#' element labels to the scene at each atomic positions. \code{info3d} activate 
#' an interactive mode to add labels by selecting atoms by \bold{right-clicking} 
#' on the current \sQuote{rgl} scene. To escape the interactive mode press the 
#' ESC key. The labels are as follow: "ResidResname:EleidElename"
#' 
#' @return \code{addResLab} and \code{addEleLab} return (using invisible) a
#' two-column data.frame containing the IDs and type indicators of the objects
#' added to the scene.
#' 
#' @param x an R object containing atomic coordinates.
#' @param at.centre a single element logical vector indicating if residue labels
#'   have to be added only at the position of the residue's centre-of-mass 
#'   instead of at each atomic position.
#' @param eleid a single element logical vector indicating if the element ids 
#'   have to be concatenated with the element names to prepare the labels.
#' @param col the colors used to display the labels.
#' @param adj one value specifying the horizontal adjustment, or two, specifying 
#'   horizontal and vertical adjustment respectively. See 
#'   \code{\link{text3d}}
#' @param id vector of ID numbers of \sQuote{rgl} items, as returned by 
#'   \code{rgl.ids}. The vertexes of these items are used to display the labels.
#' @param verbose a logical value specifying if information have to be printed 
#'   to the terminal.
#' @param \dots further arguments passed to or from other methods.
#'   
#' 
#' @seealso \code{\link{pdb}}, \code{\link{visualize}}, \code{\link{measure}}
#'
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' visualize(x, type = "l", mode = NULL)
#' addResLab(x)
#' x <- read.pdb(system.file("examples/Pentacene.pdb", package="Rpdb"))
#' visualize(x, type = "l", mode = NULL)
#' addEleLab(x)
#' 
#' @keywords dynamic
#'       
#' @name addLabels
#' @export
addResLab <- function(x, ...)
  UseMethod("addResLab")

#' @rdname addLabels
#' @export
addResLab.atoms <- function(x, at.centre = TRUE, col = "black", ...){
  if(at.centre){
    labels <- paste0(x$resname, x$resid)
    labels <- labels[!duplicated(labels)]
    xyz <- centres(x)
  }
  else{
    labels <- paste0(x$resname, x$resid)
    xyz <- coords(x)
  }
  txt.ids <- rgl::text3d(xyz, text=labels, col = col, ...)
  invisible(txt.ids)
}

#' @rdname addLabels
#' @export
addResLab.pdb <- function(x, at.centre = TRUE, col = "black", ...)
  addResLab.atoms(x$atoms, at.centre = at.centre, col = col, ...)

#' @rdname addLabels
#' @export
addEleLab <- function(x, ...)
  UseMethod("addEleLab")

#' @rdname addLabels
#' @export
addEleLab.atoms <- function(x, eleid = FALSE, col = "black", ...){
  if(eleid){
    labels <- paste(x$elename, x$eleid, sep=":")
  }
  else{
    labels <- x$elename    
  }
  xyz <- coords(x)
  txt.ids <- rgl::text3d(xyz, text=labels, col = col, ...)
  invisible(txt.ids)
}

#' @rdname addLabels
#' @export
addEleLab.pdb <- function(x, eleid = FALSE, col = "black", ...)
  addEleLab.atoms(x$atoms, eleid = eleid, col = col, ...)

#' @rdname addLabels
#' @export
info3d <- function(...)
  UseMethod("info3d")

#' @rdname addLabels
#' @export
info3d.atoms <- function(x, id = rgl::rgl.ids(), col = "black", verbose = TRUE, adj = 0, ...){
  cat("Presse ESC when you have finish your selections.\n")
  ids <- id[id$type!="text",]$id
  if(basis(x) == "abc")
    stop("'x' must contain Cartesian coordinates")
  xyz <- coords(x, basis = "xyz")
  
  sph.ids <- NULL
  txt.ids <- NULL
  clean <- function(){
    if(!is.null(sph.ids))
      rgl::rgl.pop(id = sph.ids)
    if(!is.null(txt.ids))
      rgl::rgl.pop(id = txt.ids)
  }
  on.exit(clean())
  
  dist <- 0.0015
  sel <- NULL
  while(TRUE){
    f <- rgl::rgl.select3d(button = "right", ...)
    if(is.null(f)) 
      break
    e <- environment(f)
    rgl.info <- lapply(ids,
                       function(id){
                         verts <- rgl::rgl.attrib(id,"vertices")
                         info <- verts
                         return(info)
                       })
    info <- do.call(rbind, rgl.info)
    verts <- info[,1:3]
    hits <- f(verts)
    if(!any(hits)){
      wincoords <- rgl::rgl.user2window(verts, projection = e$proj)
      hits <- (0 <= wincoords[,3]) && (wincoords[,3] <= 1)
      if(any(hits)){
        dists <- (wincoords[, 1] - e$llx)^2 +
          (wincoords[, 2] - e$lly)^2
        hits <- (dists < dist) & (dists == min(dists))
      }
    }
    if(any(hits)){
      verts <- unique(verts[hits, , drop = FALSE])
      sel <- rbind(sel, verts)
      res.lab <- do.call(paste0, x[,c("resid","resname")])
      ele.lab <- do.call(paste0, x[,c("eleid","elename")])
      labels <- paste(res.lab, ele.lab, sep=":")
      M <-     xyz$x1%in%round(verts[,1],digits=3)
      M <- M & xyz$x2%in%round(verts[,2],digits=3)
      M <- M & xyz$x3%in%round(verts[,3],digits=3)
      if(verbose){
        if(length(which(M))!=0)
          print(subset(x, M))
      }
      if(length(which(M))!=0){
        txt.ids <- c(txt.ids, 
                     rgl::text3d(verts, text=labels[M], col = col, adj = adj, ...))
      }
    }
    else{
      clean()
      sph.ids <- NULL
      txt.ids <- NULL
      sel <- NULL
    }
  }
}

#' @rdname addLabels
#' @export
info3d.pdb <- function(x, id = rgl::rgl.ids(), col = "black", verbose = TRUE, adj = 0, ...)
  info3d.atoms(x$atoms, id = id, col = col, verbose = verbose, adj = adj, ...)
