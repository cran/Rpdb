
### Various Helper Functions


### Check Basis
check.basis = function(x) {
	if( ! x %in% c("xyz","abc"))
		stop("'basis' must be equal to 'xyz' or 'abc'");
}

# id or Logical
check.idOrLogical = function(x, max, lbl) {
	if(any(round(x) != x))
		stop(paste0(lbl, " must be a logical or integer vector"));
    if(any(x > max))
      stop(paste0(lbl, " contains indices out of range"));
}

