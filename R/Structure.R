
### Molecular Structure


is.structure = function(x) {
	if(inherits(x, "structure")) return(TRUE);
	return(FALSE);
}

as.pdb.structure = function(helix = NULL, sheet = NULL, ...) {
	x = list(Helix = helix, Sheet = sheet);
	class(x) = c("structure", "list");
	return(x);
}

# Note:
# - Positions are coded using 4 digits;
# - Duplicate positions are disambiguated using the Insrtion Code;
# - Insertion codes use letters: potential for 26*2 more positions;
# - Sense of Strand:
#   0 = First strand;
#   1 = Parallel to previous, -1 = Anti-parallel to previous;
read.pdb.structure = function(x) {
	### Helix:
	helix = read.pdb.helix(x);
	### Sheet:
	sheet = read.pdb.sheet(x);
	#
	structMol = as.pdb.structure(helix, sheet);
	return(structMol);
}

###
# Helix Types:
# 1 = Right-handed Alpha;  2 = Right-handed Omega;
# ...
read.pdb.helix = function(x) {
	# https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX
	isHelix = grepl("^HELIX ", x);
	idHelix = which(isHelix);
	isHelix = length(idHelix) > 0;
	if(! isHelix) return(NULL);
	# Helix:
	sH = x[idHelix];
	hID   = trim(substr(sH,  8, 10));
	hName = trim(substr(sH, 12, 14));
	# Start:
	hAAS  = trim(substr(sH, 16, 18)); # Residue name
	hChS  = trim(substr(sH, 20, 20)); # Chain ID
	hAASn = trim(substr(sH, 22, 25)); # Seq No.
	hAASI = trim(substr(sH, 26, 26)); # Insert code
	# End:
	hAAE  = trim(substr(sH, 28, 30)); # Residue name
	hChE  = trim(substr(sH, 32, 32)); # Chain ID
	hAAEn = trim(substr(sH, 34, 37)); # Seq No.
	hAAEI = trim(substr(sH, 38, 38)); # Insert code
	hType = trim(substr(sH, 39, 40)); # Helix Class
	hComm = trim(substr(sH, 41, 70)); # Comments
	hLen  = trim(substr(sH, 72, 76)); # Length Helix
	helix = data.frame(ID = as.integer(hID), Name  = hName,
		Type = as.integer(hType), # Type of Helix (including handedness)
		Len = as.integer(hLen),
		chS = hChS, chE = hChE,
		AAS = hAAS, idAAS = as.integer(hAASn), insertS = hAASI,
		AAE = hAAE, idAAE = as.integer(hAAEn), insertE = hAAEI,
		Comment = hComm
		);
	return(helix);
}
read.pdb.sheet = function(x) {
	# https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET
	isSheet = grepl("^SHEET ", x);
	idSheet = which(isSheet);
	isSheet = length(idSheet) > 0;
	if(! isSheet) return(NULL);
	# Sheets:
	sSh = x[idSheet];
	sID   = trim(substr(sSh,  8, 10)); # Strand ID
	sName = trim(substr(sSh, 12, 14)); # Sheet ID
	sNStr = trim(substr(sSh, 15, 16)); # No. Strands
	# Start:
	sAAS  = trim(substr(sSh, 18, 20)); # AA Code
	sChS  = trim(substr(sSh, 22, 22)); # Chain ID
	sAASn = trim(substr(sSh, 23, 26)); # Seq No.
	sAASI = trim(substr(sSh, 27, 27)); # Insertion Code
	# End:
	sAAE  = trim(substr(sSh, 29, 31)); # AA Code
	sChE  = trim(substr(sSh, 32, 33)); # Chain ID
	sAAEn = trim(substr(sSh, 34, 37)); # Seq No
	sAAEI = trim(substr(sSh, 38, 38)); # Insertion Code
	sType = trim(substr(sSh, 39, 40)); # Sense of Strand
	### Registration
	# Current Strand:
	atomRc   = trim(substr(sSh, 42, 45)); # Atom
	insertRc = trim(substr(sSh, 55, 55)); # Insertion
	# Previous Strand:
	atomRp   = trim(substr(sSh, 42, 45)); # Atom
	insertRp = trim(substr(sSh, 70, 70)); # Insertion
	# TODO: remaining fields;
	sheet = data.frame(ID = as.integer(sID), Name  = sName,
		Type = as.integer(sType), # Sense of Strand
		chS = sChS, chE = sChE,
		AAS = sAAS, idAAS = as.integer(sAASn), insertS = sAASI,
		AAE = sAAE, idAAE = as.integer(sAAEn), insertE = sAAEI,
		# Registration
		atomRc   = atomRc,
		insertRc = insertRc,
		atomRp   = atomRp,
		insertRp = insertRp
		);
	return(sheet);
}
