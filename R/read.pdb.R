#' PDB File Reader
#' 
#' Reads a Protein Data Bank (PDB) coordinates file.
#' 
#' The \code{read.pdb} function reads the TITLE, REMARK, ATOM, HETATM, CRYST1 and CONECT records from a PDB file. Three different reading modes can be used depending on the value of \code{MODEL}: 
#' \itemize{
#'   \item When \code{MODEL} is a vector of integers, MODEL sections whose serial numbers match these integers are read.
#'   \item When \code{MODEL == NULL}, all MODEL sections are read.
#'   \item When \code{MODEL == NA}, MODEL records are ignored and all ATOM and/or HETATM records are merged together to return a single object.
#' }
#' When multiple models are specified, each of the models is actually stored
#'   as a separate pdb molecule in a list of pdb molecules.
#' If the \code{resolution} parameter is set, the function attempts to extract the resolution
#'   from the REMARKS field. Note: The resolution is only meaningfull for X-ray crystallography.
#' 
#' @return 
#' When a single MODEL section is read, this function returns an object of class  \sQuote{pdb} (a list with a \code{class} attribute equal to \code{pdb}) with the following components:
#' \item{title}{a character vector containing the TITLE records found in the PDB file.}
#' \item{remark}{a character vector containing the REMARK records found in the PDB file.}
#' \item{crystal}{a list of class \sQuote{crystal} containing the first CRYST1 record found in the PDB file. All others are ignored.}
#' \item{atoms}{a data.frame of class \sQuote{atoms} containing the ATOM and HETATM records found in the PDB file.}
#' \item{conect}{a data.frame of class \sQuote{conect} containing the CONECT records found in the PDB file.}
#' When multiple MODEL sections are read, a list of object of class \sQuote{pdb} is returned.
#' 
#' @param file a single element character vector containing the name of the PDB file to be read.
#' @param ATOM a logical value indicating whether to read the ATOM records.
#' @param HETATM a logical value indicating whether to read the HETATM records.
#' @param CRYSTAL a logical value indicating whether to read the crystal cell parameters (from the CRYST1 pdb record).
#' @param CONECT a logical value indicating whether to read the CONECT records.
#' @param TITLE a logical value indicating whether to read the TITLE records.
#' @param REMARK a logical value indicating whether to read the REMARK records.
#' @param MODEL an integer vector containing the serial number of the MODEL sections to be read. Can also be equal to NULL to read all the MODEL sections or to NA to ignore MODEL records (see details).
#' @param CRYST1 will be replaced by the CRYSTAL argument; existing code should be migrated to use the CRYSTAL argument.
#' @param resolution logical value indicating wheter to extract the resolution (see details).
#' @param verbose logical value indicating wheter to print additional information, e.g. number of models.
#' 
#' @references 
#' PDB format has been taken from:
#' http://www.wwpdb.org/documentation/format33/v3.3.html
#' 
#' @seealso
#' \code{\link{write.pdb}}, \code{\link{pdb}}, \code{\link{crystal}}, \code{\link{atoms}}, \code{\link{conect}}
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
                     CONECT = TRUE, TITLE = TRUE, REMARK = TRUE, MODEL = 1,
					 CRYST1 = NULL, resolution = TRUE, verbose = TRUE)
{
	if( ! is.null(CRYST1)) {
		if(CRYSTAL) {
			# Temporary:
			CRYSTAL = CRYST1;
		} else {
			stop("Please provide only a logical value for CRYSTAL!",
			"Argument CRYST1 will be deprecated!");
		}
	}
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
	if(TITLE & any(isTitle)) {
		title = subset(lines, isTitle);
		if(length(title) == 1) {
		# Text starts actually at npos = 11 as well;
		title = substr(title, 7, 80);
		} else {
			title = substr(title, 11, 80);
		}
		title = trim(title);
	}
	
	### Remarks:
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
		dfResolution = extract.pdb.Resolution(tmp);
	}
  
  ### Atoms:
  atom <- character(0)
  is.hetatm <- rep(FALSE, length(recname))
  is.atom   <- rep(FALSE, length(recname))
  if(ATOM  ) is.atom   <- recname == "ATOM  "
  if(HETATM) is.hetatm <- recname == "HETATM"
  atoms <- subset(lines, is.atom | is.hetatm)
  if(length(atoms) == 0) stop("No atoms have been selected")
  
  ### Models:
  # NUMMDL: is optional;
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
  model.factor <- rep(0, length(recname))
  model.ids   <- "MODEL.1"
  model.start <- which(recname == "MODEL "); # grep("^MODEL ", recname)
  model.end   <- which(recname == "ENDMDL"); # grep("^ENDMDL", recname)
  if(length(model.start) != length(model.end))
    stop("'Unterminated MODEL section'");
  hasModel = ! (length(model.start) == 0);
  if(hasModel) {
    if(any(model.start >= model.end))
      stop("'Unterminated MODEL section'");
	# All models:
    if(is.null(MODEL)) {
      model.ids <- as.integer(substr(lines[model.start], 11, 14))
      MODEL <- model.ids 
    }
    if( ! is.na(MODEL[1])) {
      model.ids <- as.integer(substr(lines[model.start], 11, 14));
	  if(verbose) {
        if(length(model.ids) > length(MODEL))
		  cat("Note: loading only models: ",
		    paste0(MODEL, collapse = ", "), ";\n", sep = "");
	  }
	  idModels = which(model.ids %in% MODEL);
      model.start <- model.start[idModels]
      model.end   <- model.end  [idModels]
      model.ids   <- model.ids  [idModels]
      model.ids   <- paste("MODEL", model.ids, sep=".");
      # model.factorE = model.factor; # just rep(0, ...)
      # model.factor [model.start + 1] = model.start;
      # model.factorE[model.end      ] = model.start;
      # model.factor = cumsum(model.factor) - cumsum(model.factorE);
	  # rm(model.factorE);
		model.factor [model.start + 1] = model.start;
		model.factor [model.end      ] = - model.start;
		model.factor = cumsum(model.factor);
		model.factor[model.factor == 0] <- NA
    }
  } else if(verbose && nModels > 0) {
    cat("Number of actual models: 0\n");
  }
  model.factor <- as.factor(model.factor)
  levels(model.factor) <- model.ids
  model.factor <- model.factor[is.atom | is.hetatm]
  
  ### Atoms:
  recAtom <- trim(substr(atoms,  1,  6))
  eleid   <- trim(substr(atoms,  7, 11))
  elename <- trim(substr(atoms, 13, 16))
  alt     <- trim(substr(atoms, 17, 17))
  resname <- trim(substr(atoms, 18, 21))
  chainid <- trim(substr(atoms, 22, 22))
  resid   <- trim(substr(atoms, 23, 26))
  insert  <- trim(substr(atoms, 27, 27))
  x1      <-      substr(atoms, 31, 38)
  x2      <-      substr(atoms, 39, 46)
  x3      <-      substr(atoms, 47, 54)
  occ     <-      substr(atoms, 55, 60)
  temp    <-      substr(atoms, 61, 66)
  segid   <- trim(substr(atoms, 73, 75))

  atoms <- atoms.default(recAtom, eleid, elename, alt,
                      resname, chainid, resid, insert,
                      x1, x2, x3, occ, temp, segid, basis = "xyz")
  
  ### Crystal Cell:
  cryst1 <- NULL;
  # Note: could also use recname;
  # isCrystal = grepl("^CRYST1", lines);
  isCrystal = (recname == "CRYST1");
  if(CRYSTAL && any(isCrystal)) {
    cryst1 <- subset(lines, isCrystal);
    cryst1 <- as.crystal.character(cryst1);
  } else if(CRYSTAL) {
    warning("No 'CRYSTAL' record!");
	cryst1 = NULL; # TODO
  }
  
  ### Connections:
  conect = NULL
  isConnect = grepl("^CONECT", lines);
  if(CONECT && any(isConnect)) {
    conect <- subset(lines, isConnect);
    C0 <- as.integer(substr(conect,  7, 11))
    C1 <- as.integer(substr(conect, 12, 16))
    C2 <- as.integer(substr(conect, 17, 21))
    C3 <- as.integer(substr(conect, 22, 26))
    C4 <- as.integer(substr(conect, 27, 31))
    C0 <- rep(C0, 4)
    C1 <- c(C1,C2,C3,C4)
    C0 <- subset(C0, !is.na(C1))
    C1 <- subset(C1, !is.na(C1))
    conect <- conect.default(eleid.1 = C0, eleid.2 = C1)
    tokeep <- conect$eleid.1 %in% atoms$eleid & conect$eleid.2 %in% atoms$eleid
    conect <- subset(conect, tokeep)
  }

  to.return <- pdb(atoms, cryst1, conect, remark, title,
    resolution = dfResolution);
  to.return <- split(to.return, model.factor)
  if(length(to.return) == 1) to.return <- to.return[[1]]

  return(to.return)
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

extract.pdb.Resolution = function(x) {
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
