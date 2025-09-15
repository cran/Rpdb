#' Centres-of-Geometry: Rolling Window
#' 
#' Computes the geometric centres using a rolling window on a polypeptide chain
#' 
#' \sQuote{centres.ppRoll} is a generic function to compute the geometric centres
#' using a rolling window on an object containing atomic coordinates for
#' a polypeptide chain. Due to limitations by CRAN check, the function is further
#' referred as \sQuote{centres}.
#'
#' The function may be useful to visualize the overall protein structure,
#' but in a less cluttered way. Unlike the protein backbone, the centres
#' may capture better the bulk of the protein.
#'
#'
#'   
#' @param x an R object containing atomic coordinates (with additional class \code{ppRoll}).
#' @param window size of the rolling window, specified as number of amino-acids;
#' @param chain apply only to the respective chains, by default all chains;
#' @param na.rm a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}.
#'
#' @return Return an object of class \sQuote{coords} containing the coordinates
#'   of the centres.
#' 
#' @examples
#' ### Example 1: Toll-Like Receptor 2
#' # Download pdb file for TLR2: 2Z81;
#' # Read the pdb-file:
#' # x = read.pdb("pdb2z81.ent")
#' 
#' # Compute the centres:
#' # tmp = centres.ppRoll(x)
#'
#' # Visualize the structure:
#' # visualize(x, type="p", pbc.box = FALSE)
#' # lines3d(tmp, lwd = 5, col = "red")
#'
#' # Add Another "strand":
#' # tmp = centres.ppRoll(x, window = 8)
#' # lines3d(tmp, lwd = 5, col = "#FA3296")


#' @name centres.ppRoll
#' @method centres ppRoll
#' @exportS3Method centres ppRoll
centres.ppRoll <- function(x, window, ...) {
	UseMethod("centres.ppRoll");
}

#' @rdname centres.ppRoll
#' @method centres.ppRoll default
#' @exportS3Method centres.ppRoll default
centres.ppRoll.default = function(x, window = 34, na.rm = TRUE, ...) {
	ids = unique(x$resid[x$elename == "CA"]);
	len = length(ids);
	if(len == 0) {
		xyz = data.frame(x = numeric(0), y = numeric(0), z = numeric(0));
		return(as.coords(xyz));
	}
	if(len <= window) {
		xyz = data.frame(
			x = mean(x$x1, na.rm=na.rm),
			y = mean(x$x2, na.rm=na.rm),
			z = mean(x$x3, na.rm=na.rm));
		return(as.coords(xyz));
	}
	ids = sort(ids);
	maxID = ids[len] - window;
	xyz = lapply(seq(ids[1], maxID, by = 1), function(id) {
		idw = which(ids >= id & ids <= id + window);
		idw = ids[idw];
		isW = x$resid %in% idw;
		xyz = data.frame(
			x = mean(x$x1[isW], na.rm=na.rm),
			y = mean(x$x2[isW], na.rm=na.rm),
			z = mean(x$x3[isW], na.rm=na.rm));
	});
	xyz = do.call(rbind, xyz);
	return(invisible(as.coords(xyz)));
}

#' @rdname centres.ppRoll
#' @method centres.ppRoll atoms
#' @exportS3Method centres.ppRoll atoms
centres.ppRoll.atoms = function(x, window = 34, chain = NULL, na.rm = TRUE, ...) {
	if(is.null(chain)) {
		chain = unique(x$chainid);
	}
	FUN = function(chain) {
		xyz = x[x$chainid == chain, c("x1","x2","x3", "elename", "resid"),
			drop = FALSE];
		xyz = centres.ppRoll.default(xyz, window=window, na.rm=na.rm, ...);
		if(nrow(xyz) > 0) xyz$chainid = chain;
		return(xyz);	
	}
	xyz = lapply(chain, FUN);
	xyz = do.call(rbind, xyz);
	invisible(xyz);
}


#' @rdname centres.ppRoll
#' @method centres.ppRoll pdb
#' @exportS3Method centres.ppRoll pdb
centres.ppRoll.pdb = function(x, window = 34, chain = NULL, na.rm = TRUE, ...) {
	tmp = centres.ppRoll.atoms(x$atoms, window=window, chain=chain,
		na.rm = na.rm, ...);
	if(nrow(tmp) <= 2) return(tmp);
	invisible(tmp);
}
