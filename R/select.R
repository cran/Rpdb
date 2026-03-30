#' Select atoms
#' 
#' Select groups of atoms from a PDB object.
#' 
#' \code{select} is a generic function, which enables
#'   the selection of element IDs from molecules.

#' @return 
#' \code{select} returns a vector with the element IDs.


#' @param x an R object containing atoms;
#' @param id an integer vector specifying the id of residues.
#' @param residue a character vector specifying the residues to select.
#' @param chain character vector specifying which chains to process;
#' @param \dots arguments passed to methods.


#' @name select
#' @export
select = function(x, id, ...)
	UseMethod("select");

#' @rdname select
#' @export
select.pdb = function(x, id, residue, chain, ...) {
	x   = x$atoms;
	ids = select.atoms(x, id=id, residue=residue, chain=chain, ...);
	return(ids);
}

# Note:
# - Selects ALL atoms (not only 1 atom per residue)!
#' @rdname select
#' @export
select.atoms = function(x, id, residue, chain, ...) {
	isE = rep(TRUE, natom(x));
	if(! missing(chain)) {
		isE = x$chain %in% chain;
	}
	if(! missing(id)) {
		isE = isE & x$resid %in% id;
	}
	if(! missing(residue)) {
		residue = toupper(residue);
		isE = isE & x$resname %in% residue;
	}
	ids = x$eleid[isE];
	return(ids);
}


# Water
which.water = function(x) {
	if(inherits(x, "pdb")) x = x$atoms;
	isHOH = x$Hetero & x$resname == "HOH";
	ids = x$eleid[isHOH];
	return(ids)
}
