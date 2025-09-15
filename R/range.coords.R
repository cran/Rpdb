#' Range of Atomic Coordinates
#' 
#' Determines the range of atomic coordinates.
#' 
#' 
#' @return Return a \code{\link{data.frame}} whose columns contain the range of
#'   the first, second and third coordinates of \code{x}.
#'   
#' @param x an R object containing atomic coordinates.
#' @param na.rm logical, indicating if \code{\link{NA}}'s should be omitted.
#' @param finite logical, indicating if all non-finite elements should be omitted.
#' @param \dots further arguments passed to or from other methods. 
#'
#' @seealso \code{range}, \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}
#' 
#' @examples 
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' range(x)
#' range(range(x))
#' 
#' @keywords manip
#' 
#' @name range.coords
#' @export
range.coords <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")

  to.return <- lapply(x[,c("x1","x2","x3")], range, na.rm, finite)
  to.return <- as.data.frame(to.return, row.names = c("min","max"))

  c.names <- unlist(strsplit(basis(x), ""))
  colnames(to.return) <- c.names
  
  return(to.return)
}

#' @rdname range.coords
#' @export
range.atoms <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  to.return <- range.coords(coords(x), na.rm, finite)
  return(to.return)
}

#' @rdname range.coords
#' @export
range.pdb <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  to.return <- range.atoms(x$atoms, na.rm, finite)
  return(to.return)
}
