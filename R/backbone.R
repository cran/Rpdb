


isProtein = function(x) {
	isPr = "CA" %in% x$atoms$elename;
	return(isPr);
}

# TODO: compute backbone for each chain;
asBackbone = function(x) {
	idBB = tapply(x$atoms[, c("eleid", "elename")],
				list(x$atoms$segid, x$atoms$resid),
				function(tmp) {
			id = match(c("N", "CA", "C", "O"), tmp$elename);
			id = tmp$eleid[id];
			id = c(id, id[3]);
		})
	idBB = unlist(idBB);
	idBB = data.frame(idBB[- length(idBB)], idBB[-1]);
	names(idBB) = c("eleid.1", "eleid.2");
	connect.default(idBB);
}
