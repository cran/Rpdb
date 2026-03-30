#' Extract Chain-Codes from a PDB Molecule

#' Extract the codes of all chains from a PDB object.

#' \code{chains} is a generic function to extract the ids
#' of the chains from a PDB molecule.


#' @param x Object of class \code{pdb} or class \code{atoms};
#' @param sorted logical indicating if the returned chain-codes should be sorted;
#' @param \dots arguments passed to other methods.

#' @return \code{chains} returns the codes of the individual chains in the molecule.


#' @name chains
#' @export
chains = function(x, ...)
	UseMethod("chains")

#' @rdname chains
#' @export
chains.atoms = function(x, sorted = TRUE, ...) {
	ch = unique(x$chainid);
	if(sorted) ch = sort(ch);
	return(ch);
}

#' @rdname chains
#' @export
chains.pdb = function(x, ...) {
	chains(x$atoms);
}
