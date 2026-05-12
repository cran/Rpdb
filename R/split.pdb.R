#' Divide and Reassemble \sQuote{pdb} Objects
#' 
#' \code{split} divides a \sQuote{pdb} object by groups of atoms defined by 
#' \code{f}. \code{unsplit} reverses the effect of \code{split}.
#' 
#' \code{split} produce a list of \sQuote{pdb} objects with the same 
#' \code{crystal}, \code{title} and \code{remark} components as \code{x}. Only 
#' its \code{atoms} component is split while its \code{connect} component is 
#' cleaned to keep only the meaningful connectivity for each \sQuote{pdb} object
#' of the list returned by the function. \code{unlist} produce a \sQuote{pdb} 
#' object with the same \code{crystal}, \code{title} and \code{remark} components
#' as the first element of \code{value}. The \code{atoms} and \code{connect} 
#' components of all the elements of \code{value} are combined by row.
#' 
#' @return
#' \code{split} returns a list of \sQuote{pdb} objects containing
#'   the data for the groups of atoms. The components of the list are named
#'   by the levels of \code{f} (after converting to a factor, or
#'   if already a factor and \code{drop=TRUE}, dropping unused levels).\cr
#' \code{unsplit} returns a \sQuote{pdb} object for which
#'   \code{split(x, f)} equals \code{value}.
#' 
#' @param x an object of class \sQuote{pdb} to be divided into groups.
#' @param f a \sQuote{factor} in the sense that
#'   \code{\link{as.factor}(f)} defines the grouping, or a list of such factors
#'   in which case their interaction is used for the grouping.
#' @param drop logical indicating if levels that do not occur should be
#'   dropped (if \code{f} is a \code{factor} or a list).
#' @param value a list of 'pdb' objects compatible with a splitting of
#'   \code{x}. Recycling applies if the lengths do not match.
#' @param \dots further potential arguments passed to methods.
#' 
#' @seealso \code{\link[base]{split}}, \code{\link{unsplit}}, \code{\link{pdb}}
#' 
#' @examples 
#' \donttest{
#' ### Split a pdb file by residue IDs and write them into separate files
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#'
#' file.names <- paste0(x$atoms$resname, "_", x$atoms$resid, ".pdb")
#' file.names <- unique(file.names)
#' pdb.resid  <- split(x, x$atoms$resid)
#' length(pdb.resid)
#' useless <- mapply(write.pdb, pdb.resid, file.names)
#'
#' # Cleanup
#' unlink(file.names)
#' }
#' 
#' @keywords category
#' 
#' @name split.pdb
#' @export
split.pdb <- function(x, f, drop = FALSE, ...)
{
	if(! is.pdb(x)) stop("'x' must be an object of class 'pdb'");
	
	# Atoms:
	atoms = split(x$atoms, f, drop=drop);
	# PDB:
	to.return = lapply(atoms, pdb, x$crystal, x$connect,
		title = x$title, remark = x$remark, hetero = x$Hetero,
		structure = x$Structure, resolution = x$Resolution);
	# Connect:
  to.return <- lapply(to.return, function(x) {
		r =     x$connect$eleid.1 %in% x$atoms$eleid;
		r = r & x$connect$eleid.2 %in% x$atoms$eleid;
		if(any(r)) {
			x$connect = x$connect[r,];
		} else x["connect"] = list(NULL);
		return(x)
		}
  )
  
  return(to.return)
}

#' @rdname split.pdb
#' @export
unsplit.pdb <- function(value, f, drop = FALSE, ...)
{
	if(!all(unlist(lapply(value, is.pdb))))
		stop("'value' must be a list containing only 'pdb' objects");
  
	title   = value[[1]]$title;
	remark  = value[[1]]$remark;
	hetero  = value[[1]]$Hetero;
	crystal = value[[1]]$crystal;
	structure  = value[[1]]$Structure;
	resolution = value[[1]]$Resolution;
  
	atoms = lapply(value, function(x) return(x$atoms));
	atoms = unsplit(atoms, f, drop=drop);
	
	### Connect:
	connect = lapply(value, function(x) return(x$connect));
	connect = do.call(rbind, connect);
	lenConn = nrow(connect);
	if(lenConn > 0) rownames(connect) = 1:lenConn;
	
	### PDB:
	to.return = pdb(atoms, crystal, connect,
		title = title, remark = remark, hetero = hetero,
		structure = structure, resolution = resolution);
	
	return(to.return);
}
