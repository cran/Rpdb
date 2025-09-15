#' Visualize a Molecular Structure
#' 
#' Use the rgl library to visualize in 3D a molecular structure.
#' 
#' Three different visualization styles are allowed.
#' \itemize{
#'   \item   When \code{type="p"}: Points are drawn at each atomic positions (very light visualization mode).
#'   \item   When \code{type="l"}: Lines are drawn between bonded atoms. The connectivity of the system has to be specifyed.
#'   \item   When \code{type="s"}: Spheres are drawn at each atomic positions (heavy visualization mode).
#'   
#'   The radii of the spheres are given by \code{radii}.
#'   \itemize{
#'     \item   When \code{radii="rcov"}: Covalent radii, taken from the \code{elements} data set, are used.
#'     \item   When \code{radii="rvdw"}: Van der Waals radii, taken from the \code{elements} data set, are used.
#'     \item   When \code{radii} is a numeric vector: The numeric values are used to assign to each atom a radius. If \code{length(radii) != natom(pdb)} \code{radii} is recycled.
#'   }
#' }
#' When \code{xyz}, \code{abc} or \code{pbc.box} are \code{NULL}, the axis or pbc box are are added depending if a \sQuote{crystal} object can be found.\cr
#' Two different interactive visualization modes are avalable:
#'   \itemize{
#'     \item When \code{mode="measure"}: bond lengths, angles and dihedrals can be measured by \bold{right-clicing} on the atoms.
#'     \item When \code{mode="info"}: atomic labels can be added to the scene by \bold{right-clicing} on the atoms. The labels are as follow: "ResidResname:EleidElename"
#'   }
#' When \code{mode=NULL} the interactive mode is disabled. To escape the interactive mode press the ESC key.
#' 
#' @return 
#' Return (using invisible) a two-column data.frame containing the IDs
#' and type indicators of the objects added to the scene.
#' 
#' @param x an object or the name of a PDB file containing the molecular structure to visualize.
#' @param elename a character vector containing the atomic names used to chose atom colors and radii.
#' @param cryst1 an object of class \sQuote{crystal}. See \code{\link{crystal}}
#' @param conect an object of class \sQuote{conect}. See \code{\link{conect}}
#' @param mode a single element character vector indicating the visualization mode (See details).
#' @param type a character string indicating the visualization style (See details).
#' @param xyz a logical value indicating whether the x, y and z axes have to be added to the scene. See details
#' @param abc a logical value indicating whether the a, b and c axes have to be added to the scene. See details
#' @param pbc.box a logical value indicating whether the pbc box has to be added to the scene. See details
#' @param lwd a numeric value indication the line width used to plot the axes, the pbc box and atomic bonds when \code{type = "l"} (see details).
#' @param lwd.xyz a numeric value indicating the line width used to plot the x, y and z axes.
#' @param lwd.abc a numeric value indicating the line width used to plot the a, b and c axes.
#' @param lwd.pbc.box a numeric value indicating the line width used to plot the pbc box.
#' @param cex.xyz a numeric value indicating the magnification used to plot the labels of the x, y and z axes.
#' @param cex.abc a numeric value indicating the magnification used to plot the labels of the a, b and c axes.
#' @param col a vector indicating the colors to use to plot each atom.
#' @param bg the color of the background
#' @param radii either a character string indicating the type of radii or a numeric vector specifying the radii of each atom to use to plot atoms as spheres (see details).
#' @param add a logical value indicating whether the plot has be to added to a existing scene (see \code{rgl.cur} and \code{open3d}).
#' @param windowRect a vector of four integers indicating the left, top, right and bottom of the displayed window in pixels (see \code{par3d}).
#' @param FOV the field of view. This controls the degree of parallax in the perspective view (see par3d).
#' @param userMatrix a 4 by 4 matrix describing user actions to display the scene (see \code{par3d}).
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{addXYZ}}, \code{\link{addABC}}, \code{\link{addPBCBox}}, \code{par3d}, \code{select3d}, \code{measure}, \code{info3d}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' visualize(x, type = "l", mode = NULL)
#' visualize(x, type = "s", radii = "rcov", mode = NULL)
#' visualize(x, type = "s", radii = "rvdw", mode = NULL)
#' visualize(x, type = "p", mode = NULL)
#' visualize(subset(x, resid != 1), type = "l", mode = NULL)
#' visualize(subset(x, resid == 1), type = "s", add = TRUE, mode = NULL)
#' 
#' @keywords dynamic
#' 
#' @name visualize
#' @export
visualize <- function(...)
  UseMethod("visualize")

#' @rdname visualize
#' @export
visualize.coords <- function( x, elename = NULL, cryst1 = NULL, conect = NULL, mode = NULL,
                              type = "l", xyz = NULL, abc = NULL, pbc.box = NULL,
							  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                              cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",
							  radii = "rvdw", add = FALSE, windowRect = c(0,0,800,600), FOV = 0,
							  userMatrix = diag(4), ...){
  crystal = cryst1;
  if(!is.coords(x)) stop("'x' must be an object of class coords.")
  
  if(basis(x) == "abc") x <- abc2xyz(x, crystal);

  # Unrecognized elements are considered as dummy atoms
  if(is.null(elename)){
    warning("'elename' was not specifyed. All atom have been considered as dummy atoms.")
    elename <- rep("Xx", natom(x))
  }
  
  # symb <- toSymbols(elename, na = "Xx");
  # symb[is.na(symb)] <- "Xx"
  # M <- match(symb, Rpdb::elements[, "symb"])
  M = match.element.character(elename, na = 1);
  
  if(is.null(col)){
    col <- Rpdb::elements[M, c("red","green","blue")]
    col <- do.call(rgb, col)
  }
  if(length(col) != natom(x)) col <- rep(col, length = natom(x))
  
  if(!add){
    rgl::open3d()
    rgl::par3d(windowRect = windowRect, userMatrix=userMatrix, FOV = FOV, ...)
    rgl::bg3d(color=bg)
  }
  ids <- rgl::ids3d()
  par.save <- rgl::par3d(skipRedraw=TRUE)
  
  if(is.null(xyz) && is.null(crystal))
    xyz <- TRUE
  else
    xyz <- FALSE
  
	if(is.null(abc)) {
		if(is.null(crystal))
			abc <- FALSE
		else
			abc <- TRUE
	} else if(! is.logical(abc)) {
		warning("Invalid abc-parameter!");
		abc = FALSE;
	}
	
	if(is.null(pbc.box)) {
		if(is.null(crystal))
			pbc.box <- FALSE
		else
			pbc.box <- TRUE
	} else if(! is.logical(pbc.box)) {
		warning("Invalid abc-parameter!");
		pbc.box = FALSE;
	}
  
  if(abc & is.null(crystal)) {
    warning("Cannot find periodical boundary conditions")
    abc <- FALSE
  }
  if(pbc.box & is.null(crystal)) {
    warning("Cannot find periodical boundary conditions")
    pbc.box <- FALSE
  }
  
  if(xyz) ids <- rbind(ids, addXYZ(lwd = lwd.xyz, cex = cex.xyz))
  if(abc) ids <- rbind(ids, addABC(crystal, lwd = lwd.abc, cex = cex.abc))
  if(pbc.box) ids <- rbind(ids, addPBCBox(crystal, lwd = lwd.pbc.box))
  # print(str(ids))
  
  if(type == "l") {
    if(is.null(conect)){
      warning("Undefined connectivity: Computing connectivity from coordinates...")
      conect <- conect(x)
    }
    ind <- t(conect)
    seg.id <- rgl::segments3d(
      x$x1[ind],
      x$x2[ind],
      x$x3[ind],
      color = col[ind], lwd=lwd, ...
    )
    seg.id <- data.frame(id = seg.id, type = "atom.seg")
    ids <- rbind(ids, seg.id)
  } else if(type == "p") {
    pts.id <- rgl::points3d(
      x$x1,
      x$x2,
      x$x3,
      color = col, ...)
    pts.id <- data.frame(id = pts.id, type = "atom.pts")
    ids <- rbind(ids, pts.id)
  }
  
  if(type == "s"){
    if(is.character(radii[1])){
      if(! radii[1] %in% c("rcov", "rbo", "rvdw") )
        stop("'radii' must be one of 'rcov', 'rbo', 'rvdw' or a numerical vector")
      radii <- Rpdb::elements[M, radii[1]]
    }
    if(all(radii == 0)){
      warning("All atoms are dummy atoms. 'radii' have been set to 1")
      radii <- rep(1, natom(x))
    }
    sph.id <- rgl::spheres3d(
      x$x1,
      x$x2,
      x$x3,
      color = col, radius = radii, ...)
    sph.id <- data.frame(id = sph.id, type = "atom.sph")
    ids <- rbind(ids, sph.id)
  }

  rgl::par3d(par.save)
  if(!is.null(mode)){
    if(mode == "measure"){
      measure(x)
    }
    else if(mode == "info"){
      stop("No 'info' mode for object of class 'coords'")
    }
    else{
      stop("Unrecognized visualization mode")
    }
  }
  
  invisible(ids)
}

#' @rdname visualize
#' @export
visualize.data.frame <- function( x, elename = NULL, cryst1 = NULL, conect = NULL, mode = NULL,
                                  type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                                  lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                                  cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                                 add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  if(is.null(basis(x))){
    basis(x) <- "xyz"
    warning("No basis attribute were found. Coordinates are assumed to Cartesian.")
  }
  
  visualize(coords(x), elename, cryst1, conect, mode, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg,  radii, add, windowRect, FOV, userMatrix, ...)
}

#' @rdname visualize
#' @export
visualize.matrix <- function( x, elename = NULL, cryst1 = NULL, conect = NULL, mode = NULL,
                              type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                              lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                              cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                              add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  if(is.null(basis(x))){
    basis(x) <- "xyz"
    warning("No basis attribute were found. Coordinates are assumed to Cartesian.")
  }
  
  visualize(coords(x), elename, cryst1, conect, mode, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg,  radii, add, windowRect, FOV, userMatrix, ...)
}

#' @rdname visualize
#' @export
visualize.atoms <- function( x, cryst1 = NULL, conect = NULL, mode = NULL,
                             type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                             lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                             cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                             add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  
  ids <- visualize(coords(x), x$elename, cryst1, conect, mode=NULL, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg, radii, add, windowRect, FOV, userMatrix, ...)
  
  if(!is.null(mode)){
    if(mode == "measure"){
      measure(x)
    }
    else if(mode == "info"){
      info3d(x)
    }
    else{
      stop("Unrecognized visualization mode")
    }
  }
  
  invisible(ids)
}

#' @rdname visualize
#' @export
visualize.pdb <- function(x, mode = NULL, type = "l", xyz = NULL, abc = NULL, pbc.box = NULL,
				lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
				cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2", radii = "rvdw",
				add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...) {
  
  # calls visualize.atoms;
  ids <- visualize(x$atoms, cryst1 = x$crystal, conect = x$conect, mode=NULL, type,
            xyz, abc, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
            cex.xyz, cex.abc, col, bg, radii, add, windowRect, FOV, userMatrix, ...)
  
  if(!is.null(mode)){
    if(mode == "measure"){
      measure(x)
    }
    else if(mode == "info"){
      info3d(x)
    }
    else{
      stop("Unrecognized visualization mode")
    }
  }
  
  invisible(ids)
}

#' @rdname visualize
#' @export
visualize.character <- function(x, mode = NULL, type = "l", xyz = NULL, abc = NULL, pbc.box = NULL, lwd = 2,
                                lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
                                cex.xyz = 2, cex.abc = 2, col = NULL, bg = "#FAFAD2",  radii = "rvdw",
                                add = FALSE, windowRect = c(0,0,800,600), FOV = 0, userMatrix=diag(4), ...){
  x <- read.pdb(x)
  visualize.pdb(x, mode = mode, type = type, xyz = xyz, abc = abc, pbc.box = pbc.box, lwd = lwd,
           lwd.xyz = lwd.xyz, lwd.abc = lwd.abc, lwd.pbc.box = lwd.pbc.box,
           cex.xyz = cex.xyz, cex.abc = cex.abc, col = col, bg = bg,  radii = radii,
           add = add, windowRect = windowRect, FOV = FOV, userMatrix = userMatrix, ...)
  
}
