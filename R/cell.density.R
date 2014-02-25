cell.density <- function(...)
  UseMethod("cell.density")

cell.density.default <- function(masses, volume, ...) {
  Na <- universalConstants["Na","Value"]
  d <- sum(masses)/(volume*1E-24*Na)
#   attr(d, "unit") <- "g.cm-3"
  return(d)
}

cell.density.pdb <- function(x, ...) {
  M <- masses(x)
  V <- cell.volume(x)
  d <- cell.density(M, V)
  return(d)
}
