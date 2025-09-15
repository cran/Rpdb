#' Format "Title" Field of PDB File
#' 
#' Formats the title field of a Protein Data Bank (PDB) file.
#' 
#' This function is mostly used internally, e.g. by \code{\link{write.pdb}},
#' but may be useful in other instances as well.
#'
#' 
#' @param x a vector of strings;
#' @param \dots currently not used;
#' 
#' @return properly formatted text string.
#'
#' @seealso \code{\link{write.pdb}}
#'
#' @examples
#' format.pdb.title(c("Molecule 1", "is just an example"))

#' @name format.pdb.title
#' @method format.pdb title
#' @export format.pdb.title
format.pdb.title = function(x, ...) {
	UseMethod("format.pdb.title");
}

#' @rdname format.pdb.title
#' @method format.pdb.title character
#' @exportS3Method format.pdb.title character
format.pdb.title.character = function(x, ...) {
	if(is.null(x)) {
		# Field is mandatory!
		x = format80("TITLE    NULL");
		return(x);
	}
	noHeader = substr(x, 1, 6) != "TITLE ";
	idHeader = which(noHeader);
	if(length(idHeader) > 0) {
		if(idHeader[1] == 1) {
			x[1] = format80(x[1], name = "TITLE    ");
			idHeader = idHeader[-1];
		}
		len = length(idHeader);
		if(len > 0) {
			id = as.character(idHeader);
			sp = sapply(nchar(id), function(len) {
				paste0(rep(" ", 4 - len), collapse = "");
			});
			x[idHeader] = paste0("TITLE", sp, id, " ", x[idHeader]);
			x[idHeader] = format80(x[idHeader]);
		}
	}
	return(x);
}

format80 = function(x, name = "") {
	len = 80 - nchar(x) - nchar(name);
	# TODO: Warn if len > 80?
	sp  = sapply(len, function(len) {
			if(len <= 0) "" else paste0(rep(" ", len), collapse = "");
		});
	paste0(name, x, sp);
}
