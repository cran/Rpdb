#' Add Axes or PBC Box to the \sQuote{rgl} Scene
#' 
#' Add lattice vectors, Cartesian axes or PBC box to the current \sQuote{rgl} scene.
#' 
#' \code{addABC}: Add the lattice vectors a, b and c to the current rgl device.\cr
#' \code{addXYZ}: Add the Cartesian axes x, y and z to the current rgl device.\cr
#' \code{addPBCBox}: Add a box representing the Periodic Boundary Conditions
#'   of a molecular system.\cr
#' \code{addBBox}: Add a bounding box around the PDB molecule
#'   constructed using the lattice vectors.\cr
#' \code{bbox.pdb}: Computes the bounding box around the PDB molecule.\cr
#' 
#' @return Returns invisibly a two-column data.frame containing the IDs 
#'   and type indicators of the objects added to the scene.
#'
#' \code{bbox.pdb}: Returns the coordinates of the bounding box.\cr
#'
#'
#' @param x an object of class \sQuote{crystal} containing the unit cell parameters.
#' @param lwd a numeric value indicating the line width used to draw the axes or
#'   the PBC box.
#' @param labels a logical value indicating whether to draw the labels of the axes.
#' @param scale a scalar value, or a numeric vector of length 2 or 3 or 6, or a 3x3 matrix
#'   used to scale the PBC box or the axes; a length 3 vector is converted to a diagonal matrix,
#'   while vectors of length 2 or 4 or 6 are converted to a diagonal matrix for scaling
#'   and a relative shift;
#' @param cex a numeric value indicating the magnification used to draw the 
#'   labels of the axes.
#' @param col colour used for the PBC box;
#' @param alpha a numeric value specifying the transparency of the lines;
#' @param dx a numeric value specifying the absolute length of the xyz axes-vectors;
#' @param bidir logical indicating to draw the XYZ-axes in both directions;
#' @param adj.lab numeric value performing relative adjustment of label position;
#' @param \dots further arguments passed to the \code{rgl} methods;
#'   
#' @seealso \code{\link{visualize}}, \code{\link[rgl]{rgl.open}}, \code{\link[rgl]{par3d}},
#'   \code{\link{addLabels}}
#'   
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' visualize(x, type = "l", xyz = FALSE, abc = FALSE, pbc.box = FALSE, mode = NULL)
#' addXYZ()
#' addABC(x)
#' addPBCBox(x)
#' 
#' @keywords dynamic


#' @name addAxes
#' @export
addABC <- function(x, lwd = 2, scale = NULL, labels = TRUE,
		cex = 2, adj.lab = 4, ...) {
	if(missing(x)) stop("Please specify a 'crystal' object");
	if(is.pdb(x)) x = crystal(x);
	if(! is.crystal(x)) stop("'x' must be an object of class 'crystal");
	
	cell = cell.coords(x);
	# Scale Axes:
	mAxes = rbind(
		c(0,0,0), cell[,1],
		c(0,0,0), cell[,2],
		c(0,0,0), cell[,3]);
	if(! is.null(scale)) {
		mAxes = scaleBox(scale, mAxes, cell=cell);
	}
	seg.id = rgl::segments3d(
		mAxes,
		col = c("red","red","green","green","blue","blue"),
		lwd = lwd, ...);
	seg.id = data.frame(id = seg.id, type = "abc.seg");
	### Labels:
	lab.id = NULL;
	# Behaves sub-optimally with non-orthogonal axes;
	if(labels) {
		mLab = mAxes[c(2,4,6),];
		mOri = mAxes[c(1,3,5),];
		mV   = mLab - mOri;
		an = mV[1,] / sqrt(sum(mV[1,]^2));
		bn = mV[2,] / sqrt(sum(mV[2,]^2));
		cn = mV[3,] / sqrt(sum(mV[3,]^2));
		mLab = mLab + adj.lab * cbind(an, bn, cn);
		lab.id = rgl::text3d(
			mLab,
			texts = c("a","b","c"),
			col   = c("red","green","blue"),
			cex   = cex, adj = c(0.25, 0.25));
		lab.id = data.frame(id = lab.id, type = "abc.lab");
	}
	# All:
	ids = rbind(seg.id, lab.id);
	invisible(ids);
}

#' @rdname addAxes
#' @export
addXYZ <- function(lwd = 2, scale = NULL, labels = TRUE,
		cex = 2, col = "black", alpha = NULL, ..., dx = 5, bidir = FALSE) {
	if(length(dx) == 1) dx = c(dx,dx,dx);
	if(length(dx) != 3) stop("Invalid value for dx!");
	dy = dx[2]; dz = dx[3]; dx = dx[1];
	# bi-directional: (-dx, dx);
	mAxes = rbind(
		c(0,0,0), c( dx,  0,  0),
		c(0,0,0), c(  0, dy,  0),
		c(0,0,0), c(  0,  0, dz));
	if(bidir) mAxes = rbind(mAxes, - mAxes);
	# Scale Axes:
	if(! is.null(scale)) {
		# Note: What cell value to use: = 1 or dx?
		mAxes = scaleBox(scale, mAxes, cell = diag(1, 3));
	}
	seg.id = rgl::segments3d(
		mAxes, lwd = lwd, col = col, alpha = alpha, ...);
	seg.id = data.frame(id = seg.id, type = "xyz.seg");
	# Labels:
	lab.id = NULL;
	if(labels) {
		dx = mAxes[c(2,4,6),]; # scaling is included;
		lab.id = addLabelsXYZ(cex=cex, col=col, dx=dx);
	}
	#
	ids = rbind(seg.id, lab.id);
	invisible(ids)
}

addLabelsXYZ = function(labels = c("x","y","z"), cex = 2, col = "black",
		dx = c(5,5,5), adj.lab = 0.2) {
	if(is.matrix(dx)) {
		mLab = dx %*% diag(1 + adj.lab, 3);
	} else {
		if(length(dx) == 1) dx = rep(dx, 3);
		dy = dx[2]; dz = dx[3]; dx = dx[1];
		dd = 1.0 + adj.lab;
		dx = dx + dd; dy = dy + dd; dz = dz + dd;
		mLab = rbind(
			c(0,0,0) + c( dx, 0.0, 0.0),
			c(0,0,0) + c(0.0,  dy, 0.0),
			c(0,0,0) + c(0.0, 0.0,  dz));
	}
	lab.id = rgl::text3d(
		mLab,
		texts = labels,
		cex   = cex, col = col
	)
	lab.id = data.frame(id = lab.id, type = "xyz.lab");
}

#' @rdname addAxes
#' @export
addPBCBox = function(x, scale = NULL,
		lwd = 2, col = "black", alpha = 0.25) {
	if(missing(x)) stop("Please specify a 'crystal' object");
	x = get.PDBCrystal(x);
	# Bounding Box:
	cell = cell.coords(x);
	mBox = box.coords(cell);
	# Scale Box:
	if(! is.null(scale)) {
		mBox = scaleBox(scale, mBox, cell=cell);
	}
	# PBC Box:
	seg.id = rgl::segments3d(
		mBox,
		lwd = lwd,
		col = col,
		alpha = alpha
	);
	#
	seg.id = data.frame(id = seg.id, type = "pbc.box");
	invisible(seg.id)
}


#' @rdname addAxes
#' @export
addBBox = function(x, lwd = 2, col = "black", alpha = 0.25) {
	mBox = bbox.pdb(x);
	# BBox:
	seg.id = rgl::segments3d(
		mBox,
		lwd = lwd,
		col = col,
		alpha = alpha
	);
	#
	seg.id = data.frame(id = seg.id, type = "pbc.box");
	invisible(seg.id)
}
#' @rdname addAxes
#' @export
bbox.pdb = function(x) {
	xyzB = range.lattice.pdb(x);
	xyzB = t(xyzB);
	dBox = xyzB[,2] - xyzB[,1];
	cell = cell.coords(x$crystal);
	origin = as.vector(cell %*% xyzB[, 1, drop = FALSE]);
	cell = cell %*% diag(dBox);
	bb = box.coords(cell, origin = origin);
	return(bb);
}


### Helper:

scaleBox = function(scale, mBox, cell) {
	len = length(scale);
	if(len == 1) {
		mBox = scale * mBox;
	} else if(len == 2 || len == 4 || len == 6) {
		# Diag & Relative Shift:
		if(len == 6) {
			m1 = diag(scale[1:3]);
			m2 = scale[4:6];
		} else if(len == 2) {
			m1 = diag(scale[1], 3);
			m2 = rep(scale[2], 3);
		} else {
			m1 = diag(scale[1], 3);
			m2 = scale[2:4];
		}
		mBox  = mBox %*% m1;
		shift = cell %*% m2;
		shift = rep(shift, each = nrow(mBox));
		mBox = mBox + shift;
	} else {
		if(len == 3) {
			scale = diag(scale);
		} else if(len == 9) {
			# OK
		} else stop("Invalid scale argument!");
		mBox = mBox %*% scale;
	}
	return(mBox);
}

box.coords = function(x, origin = c(0,0,0)) {
	cell_12 = x[,1] + x[,2] + origin;
	cell_13 = x[,1] + x[,3] + origin;
	cell_23 = x[,2] + x[,3] + origin;
	cell_Sm = cell_12  + x[,3]; # Opposite Point
	# Shift Cell:
	x = x + rep(origin, 3);
	mBox = rbind(
		origin,   x[,1],
		origin,   x[,2],
		origin,   x[,3],
		cell_12,  x[,1],
		cell_12,  x[,2],
		cell_12,  cell_Sm,
		cell_13,  x[,1],
		cell_13,  cell_Sm,
		cell_13,  x[,3],
		cell_23,  cell_Sm,
		cell_23,  x[,2],
		cell_23,  x[,3]
	);
	rownames(mBox) = NULL;
	return(mBox);
}
