#' PDB File Reader
#' 
#' Reads a Protein Data Bank (PDB) coordinates file.
#' 
#' The \code{read.pdb} function reads the most important records from the PDB file,
#'   including TITLE, REMARK, ATOM, HETATM, Heterogen Names (HETNAM),
#'   CRYSTAL, CONNECT, HELIX and SHEET records.
#'
#' Three different reading modes can be used depending on the value of \code{MODEL}: 
#' \itemize{
#'   \item When \code{MODEL} is a vector of integers, only the models which
#'      match these integers are read.
#'   \item When \code{MODEL == NULL}, all models are read.
#'   \item When \code{MODEL == NA}, MODEL records are ignored and all ATOM and/or HETATM records
#'     are merged together to return a single object.
#' }
#' When multiple models are specified, each of the models is actually stored
#'   as a separate pdb molecule in a list of pdb molecules.
#'   Note: The handling of models may change in a future version of Rpdb.
#' \cr
#' If the \code{resolution} parameter is set, the function attempts to extract the resolution
#'   from the REMARKS field. Note: The resolution is only meaningful for X-ray crystallography.
#' 
#' @return 
#' When a single MODEL section is read, this function returns an object of class \sQuote{pdb} (a list with a \code{class} attribute equal to \code{pdb}) with the following components:
#' \item{title}{a character vector containing the TITLE records found in the PDB file.}
#' \item{remark}{a character vector containing the REMARK records found in the PDB file.}
#' \item{crystal}{a list of class \sQuote{crystal} containing the first CRYSTAL record found in the PDB file. All others are ignored.}
#' \item{atoms}{a data.frame of class \sQuote{atoms} containing the ATOM and HETATM records found in the PDB file.}
#' \item{connect}{a data.frame of class \sQuote{connect} containing the CONNECT records found in the PDB file.}
#' \item{Structure}{a list with components\code{Helix} and \code{Sheet}}
#' \item{Resolution}{a numeric value with the resolution.}
#'
#' When multiple MODEL sections are read, a list of object of class \sQuote{pdb} is returned.
#' Note: this may change in the future.
#' 
#' @param file a single element character vector containing the name of the PDB file to be read.
#' @param ATOM a logical value indicating whether to read the ATOM records.
#' @param HETATM a logical value indicating whether to read the HETATM records.
#' @param CRYSTAL a logical value indicating whether to read the crystal cell parameters (from the CRYST1 pdb record).
#' @param CONNECT a logical value indicating whether to read the CONNECT records.
#' @param TITLE a logical value indicating whether to read the TITLE records.
#' @param REMARK a logical value indicating whether to read the REMARK records.
#' @param MODEL an integer vector containing the serial number of the MODEL sections to be read. Can also be equal to NULL to read all the MODEL sections or to NA to ignore MODEL records (see details).
#' @param resolution logical value indicating whether to extract the resolution (see details).
#' @param verbose logical value indicating whether to print additional information, e.g. number of models.
#' 
#' @references 
#' PDB format has been taken from:
#' http://www.wwpdb.org/documentation/format33/v3.3.html
#' 
#' @seealso
#' \code{\link{write.pdb}}, \code{\link{pdb}}, \code{\link{crystal}},
#' \code{\link{atoms}}, \code{\link{connect}}
#' 
#' @examples 
#' # Read a PDB file included with the package
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' 
#' # Visualize the PDB file
#' visualize(x, mode = NULL)
#' 
#' # Write the 'pdb' object 'x' in file "Rpdb.pdb" into the current directory
#' write.pdb(x, file = "Rpdb.pdb")
#'
#' # Cleanup
#' unlink("Rpdb.pdb")
#' 
#' @keywords IO
#' 
#' @name read.pdb
#' @export
read.pdb <- function(file, ATOM = TRUE, HETATM = TRUE, CRYSTAL = TRUE,
                     CONNECT = TRUE, TITLE = TRUE, REMARK = TRUE, MODEL = 1,
					 resolution = TRUE, verbose = TRUE)
{
	if(!file.exists(file))
		stop("File '", file, "' is missing!");
	
	lines = readLines(file);
	recname = substr(lines, 1, 6);
	
	trim = function(str)
		sub(' +$', '', sub('^ +', '', str));
	
	if(verbose) {
		### Obsolete:
		isObsolete = (recname == "OBSLTE");
		idObsolete = which(isObsolete);
		if(length(idObsolete) > 0) {
			warning("Obsolete pdb!");
			txtObsolete = lines[idObsolete];
			txt = trim(substr(txtObsolete, 22, 80));
			txt[1] = paste0("PDB ",
				substr(txt[1],  1,  4), " replaced by: ",
				substr(txt[1], 11, 59));
			cat(txt, sep = "\n");
			rm(txtObsolete);
		}
		### Caveat:
		isCaveat = (recname == "CAVEAT");
		idCaveat = which(isCaveat);
		if(length(idCaveat) > 0) {
			warning("The pdb contains a Caveat record!");
			txtCaveat = lines[idCaveat];
			txt = trim(substr(txtCaveat, 20, 80));
			# TODO
		}
	}
	
	### Title:
	title = NULL;
	isTitle = (recname == "TITLE ");
	if(TITLE) {
		if(any(isTitle)) {
			title = subset(lines, isTitle);
			isTitle = TRUE;
		} else {
			# Only HEADER
			isTitle = (recname == "HEADER");
			title   = subset(lines, isTitle);
			isTitle = FALSE; # actually NOT a Title;
		}
		nposE = if(isTitle) 80 else 59;
		title = substr(title, 11, nposE);
		title = trim(title);
	}
	
	### Remarks:
	# TODO: improved extraction;
	remark = NULL;
	isRemark  = (recname == "REMARK");
	hasRemark = any(isRemark);
	if(REMARK && hasRemark)
		remark = subset(lines, isRemark);
	### Resolution:
	dfResolution = NULL;
	if(resolution && hasRemark) {
		tmp = if(is.null(remark)) subset(lines, isRemark) else remark;
		# can start both at +11 & + 12;
		dfResolution = read.pdb.Resolution(tmp);
	}
  
	### Atoms:
	is.hetatm = rep(FALSE, length(recname));
	is.atom   = rep(FALSE, length(recname));
	if(ATOM  ) is.atom   = recname == "ATOM  ";
	if(HETATM) is.hetatm = recname == "HETATM";
	isAtomField = is.atom | is.hetatm;
	atoms = subset(lines, isAtomField);
	if(length(atoms) == 0) stop("No atoms have been selected");
	### Atoms:
	atoms = as.atoms.character(atoms, isHetero = is.hetatm[isAtomField]);
  
	### Models:
	# Field NUMMDL: is optional;
	isNModels = (recname == "NUMMDL");
	if(any(isNModels)) {
		txtModels = subset(lines, isNModels);
		if(length(txtModels) > 1) {
			warning("Malformed pdb file: incorrect number of models!");
			txtModels = txtModels[1];
		}
		nModels = as.integer(substr(txtModels, 7, 80)); # TODO: 7, 11?
	} else nModels = 0;
	if(verbose) cat("Number of models: ", nModels, "\n");
	#
	model.factor = rep(0, length(recname))
	model.ids    = "MODEL.1"
	model.start  = which(recname == "MODEL ");
	model.end    = which(recname == "ENDMDL");
	if(length(model.start) != length(model.end)) {
		stop("'Unterminated MODEL section'");
		# TODO: consider last ATOM field;
	}
	#
	nModels0 = 0;
	hasModel = ! (length(model.start) == 0);
	if(hasModel) {
		if(any(model.start >= model.end))
			stop("'Unterminated MODEL section'");
		# All models:
		model.ids = as.integer(substr(lines[model.start], 11, 14));
		nModels0  = length(model.ids); # Total n of models;
		if(is.null(MODEL)) {
			MODEL = model.ids;
		}
		# Only requested models:
		MODEL = MODEL[! is.na(MODEL)];
		if(length(MODEL) == 0) {
			cat("No models specified!\n");
			# TODO
		} else {
			if(verbose) {
				if(length(model.ids) > length(MODEL))
					cat("Note: Loading only models: ",
						paste0(MODEL, collapse = ", "), ";\n", sep = "");
			}
			idModels = match(MODEL, model.ids);
			naModels = is.na(idModels);
			if(any(naModels)) {
				if(verbose) {
					cat("Note: Some models not found: ",
						paste(idModels[naModels], collapse = ", "), ";\n");
				}
				idModels = idModels[! naModels];
			}
			model.start = model.start[idModels];
			model.end   = model.end  [idModels];
			model.ids   = model.ids  [idModels];
			model.ids   = paste("MODEL", model.ids, sep=".");
			model.factor [model.start + 1] = model.start;
			model.factor [model.end      ] = - model.start;
			model.factor = cumsum(model.factor);
			model.factor[model.factor == 0] = NA;
		}
	} else if(verbose && nModels > 0) {
		cat("Number of actual models: 0\n");
	}
	model.factor = as.factor(model.factor);
	levels(model.factor) = model.ids;
	model.factor = model.factor[is.atom | is.hetatm];
  
	### Crystal Cell:
	crystal = NULL;
	isCrystal = (recname == "CRYST1");
	if(CRYSTAL && any(isCrystal)) {
		crystal = subset(lines, isCrystal);
		crystal = as.crystal.character(crystal);
	} else if(CRYSTAL) {
		warning("No 'CRYSTAL' record!");
		crystal = NULL; # TODO
	}
  
	### Connections:
	connect = NULL
	isConnect = grepl("^CONECT", lines);
	if(CONNECT && any(isConnect)) {
		connect = subset(lines, isConnect);
		connect = connect.character(connect, atoms);
		# Re-index: because of TER-record;
		idA = atoms$eleid;
		atoms$eleid = seq_along(idA);
		connect$eleid.1 = match(connect$eleid.1, idA);
		connect$eleid.2 = match(connect$eleid.2, idA);
	} else {
		# match(atoms$eleid, atoms$eleid);
		if(nModels == 1) {
			atoms$eleid = seq_along(atoms$eleid);
		} else {
			# TODO: cleanup/refactor;
			factM = model.factor[! is.na(model.factor)];
			for(idM in levels(factM)) {
				atoms$eleid[factM == idM] =
					seq_along(atoms$eleid[factM == idM]);
			}
		}
	}
	### Protein Structure:
	structMol = read.pdb.structure(lines);
	
	### Hetero-Molecules:
	heteroMol = read.pdb.HeteroMol(lines, recname);
	
	### PDB Object:
	pdbObj = pdb(atoms, crystal, connect,
		title = title, remark = remark, hetero = heteroMol,
		structure = structMol, resolution = dfResolution);
	if(nModels0 > 1) {
		pdbObj = split(pdbObj, model.factor);
		if(length(pdbObj) == 1) pdbObj = pdbObj[[1]];
	}
	
	return(pdbObj)
}

### Helper

trim = function(x)
	sub(' +$', '', sub('^ +', '', x));

extract.regex = function(x, pattern, gr=0, perl=TRUE, simplify=TRUE, verbose=FALSE) {
	if(inherits(x, "data.frame")) stop("x should be an array!");
	r = regexec(pattern, x, perl=perl);
	if(verbose) cat("Finished Regex.\nStarting extraction.\n");
	gr  = gr + 1;
	len = length(gr);
	if(len > 1) {
		s = lapply(seq(length(x)), function(id) {
			tmp = r[[id]];
			if(tmp[1] == -1) return(rep("", len));
			LEN = attr(tmp, "match.length");
			sapply(seq(len), function(idGr) {
				gr = gr[idGr];
				nStart = tmp[gr];
				substr(x[id], nStart, nStart - 1 + LEN[gr]);
			})
		});
		if(simplify) s = do.call("rbind", s);
	} else {
		s = lapply(seq(length(x)), function(id) {
			tmp = r[[id]];
			if(tmp[1] == -1) return("");
			nStart = tmp[gr];
			substr(x[id], nStart, nStart - 1 + attr(tmp, "match.length")[gr]);
		});
		if(simplify) s = unlist(s);
	}
	return(s)
}

read.pdb.Resolution = function(x) {
	# can start both at +11 & + 12;
	strRes = trim(substr(x, 11, 22));
	isRes  = (strRes == "RESOLUTION.");
	if(any(isRes)) {
			x = x[isRes];
			x = trim(substr(x, 23, 80));
			res  = extract.regex(x, "([0-9.]+) ([A-Z]+)", gr = c(1,2));
			unit = res[,2];
			res  = as.numeric(res[,1]);
			isNA = is.na(res);
			unit[isNA] = x[isNA];
			attr(res, "Unit") = unit;
	} else res = NULL;
	return(res);
}

read.pdb.HeteroMol = function(x, recname) {
	idHeteroMol = which(recname == "HETNAM");
	if(length(idHeteroMol) == 0) return(NULL);
	hetero  = x[idHeteroMol];
	abbrHet = trim(substr(hetero, 11, 14));
	nmsHet  = trim(substr(hetero, 15, 70));
	# TODO: nposEND == 70?
	heteroMol = data.frame(Abbr = abbrHet, Name = nmsHet);
	return(heteroMol);
}
