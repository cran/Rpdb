


# TODO: better strategy;
isProtein = function(x) {
	isPr = "CA" %in% x$atoms$elename;
	return(isPr);
}
isProteinA = function(x) {
	isPr = "CA" %in% x$elename;
	return(isPr);
}
# TODO
isNucleicAcid.atoms = function(x) {
	nc = unique(x$resname);
	isNc = any(nc %in% c("A","C","G", "T", "U"));
	return(isNc);
}
isNucleicAcid.character = function(x) {
	nc = unique(x);
	isNc = any(nc %in% c("A","C","G", "T", "U"));
	return(isNc);
}
# Special AA:
isAASpecial = function(x) {
	x %in% c("HZP", "PTR", "TYS");
}

### Backbone:
asBackbone = function(x) {
	# TODO: check continuity
	FUN = function(tmp) {
		id = match(c("N", "CA", "C", "O"), tmp$elename);
		id = tmp$eleid[id];
		id = c(id, NA, id[3]);
	}
	FUN.Nc = backbone.nucleic();
	if(inherits(x, "atoms")) {
		atoms = x;
	} else if(inherits(x, "pdb")) atoms = x$atoms;
	# Process each chain:
	ch = chains(x);
	idBB = lapply(ch, function(ch) {
		# Modified AA: may require better checking if true AA;
		isAA = atoms$Hetero == FALSE;
		isAS = isAASpecial(atoms$resname[atoms$Hetero]);
		isAA[atoms$Hetero] = isAS;
		isCh = atoms$chainid == ch & isAA;
		tmp  = atoms[isCh, c("resid", "eleid", "elename")];
		if(nrow(tmp) == 0) return(NULL);
		#
		if(isNucleicAcid.character(atoms$resname[isCh])) {
			FUN = FUN.Nc;
		}
		idBB = tapply(tmp, tmp$resid, FUN);
		idBB = unlist(idBB);
		idBB = data.frame(idBB[- length(idBB)], idBB[-1]);
		isNA = is.na(idBB[,1]) | is.na(idBB[,2]);
		idBB = idBB[! isNA, , drop = FALSE];
		names(idBB) = c("eleid.1", "eleid.2");
		return(idBB);
	});
	idBB = do.call(rbind, idBB);
	connect.default(idBB);
}

residues = function(x, chain = NULL, hetero = NULL) {
	if(inherits(x, "pdb")) x = x$atoms;
	ch = chains(x);
	if(! is.null(chains)) {
		idCh = match(chain, ch);
		isNA = is.na(idCh);
		if(any(isNA)) {
			warning("The chains: ",
				paste(chain[isNA], collapse = ", "), " do NOT exist!");
			idCh = idCh[! isNA];
		}
		ch = ch[idCh];
		if(length(ch) == 0) return();
	}
	doHetero = ! is.null(hetero);
	tmp = lapply(ch, function(ch) {
		isCh = x$chainid == ch;
		if(doHetero) isCh = isCh & (x$Hetero == hetero);
		tmp = x[isCh, c("resname", "resid"), drop = FALSE];
		tmp = unique(tmp);
		if(nrow(tmp) > 0) tmp$chainid = ch;
		return(tmp);
	})
	tmp = do.call(rbind, tmp);
	return(tmp);
}

# 5' O-CH2-C4-C3-O-PO2-O-> 3'
backbone.nucleic = function(x) {
	FUN = function(x) {
		atom = x$elename;
		atom = sub("\\*", "'", atom);
		id = match(
			c("O5'", "C5'", "C4'", "C3'", "O3'",
			"P", "O1P", "O2P",  "C2'", "C1'", "O4'"), atom);
		id = x$eleid[id];
		id = c(id[6:7], NA, id[6], id[8], NA, # PO4
			id[6], id[-(6:11)], NA,
			id[4], id[9:11], id[3], NA, # Pentose-Ring
			id[5]);
		return(id);
	}
}