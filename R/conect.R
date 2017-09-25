#' Create \sQuote{conect} Object
#' 
#' Creates an object of class \sQuote{conect} containing the IDs of bonded atoms
#' defining the connectivity of a molecular system.
#' 
#' \code{conect} is a generic function to create objects of class 
#' \sQuote{conect}. The purpose of this class is to store CONECT records from 
#' PDB files, indicating the connectivity of a molecular system.\cr The default 
#' method creates a \code{conect} object from its different components, i.e.: 
#' \code{eleid.1} and \code{eleid.2}. Both arguments have to be specified.\cr 
#' The S3 method for object of class \sQuote{coords} determine the connectivity 
#' from atomic coordinates. A distance matrix is computed, then, for each pair 
#' of atom the distance is compared to a bounding distance computed from atomic 
#' radii. If this distance is lower than the bounding distance then the atoms 
#' are assumed to be connected.\cr The S3 method for object of class 
#' \sQuote{pdb} first use element names to search for atomic radii in the 
#' \code{elements} data set. Then atomic coordinates and radii are passed to 
#' \code{conect.coords}.\cr If \code{by.block == TRUE}, a grid is defined to 
#' determined the connectivity by block. The method is slow but allow to deal 
#' with very large systems. \cr \code{is.conect} tests if an object of class 
#' \sQuote{conect}, i.e. if it has a \dQuote{class} attribute equal to 
#' \code{conect}.
#' 
#' @return \code{conect} returns a two-column data.frame of class
#' \sQuote{conect} whose rows contain the IDs of bonded atoms. The columns of
#' this data.frame are described below: \item{eleid.1}{a integer vector
#' containing the elements IDs defining the connectivity of the system.} 
#' \item{eleid.2}{a integer vector containing the elements IDs defining the
#' connectivity of the system.} \cr\cr \code{is.conect} returns TRUE if \code{x}
#' is an object of class \sQuote{coords} and FALSE otherwise.
#' 
#' @param \dots arguments passed to methods.
#' @param eleid.1 a integer vector containing the IDs of bonded atoms.
#' @param eleid.2 a integer vector containing the IDs of bonded atoms.
#' @param x an R object containing atomic coordinates.
#' @param radii a numeric vector containing atomic radii used to find neigbours.
#' @param safety a numeric value used to extend the atomic radii.
#' @param by.block a logical value indicating whether the connectivity has to be
#'   determine by block (see details).
#'   
#' @seealso \code{\link{pdb}}
#' 
#' @examples 
#' # If atom 1 is connected to atom 2, 3, 4 and 5
#' # then we can prepare the following 'conect' object:
#' x <- conect(rep(1,4),2:5)
#' print(x)
#' is.conect(x)
#' 
#' # Compute conectivity from coordinates
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb",package="Rpdb"), CONECT = FALSE)
#' x$conect
#' x$conect <- conect(x)
#' x$conect
#' 
#' @keywords classes
#' 
#' @name conect
#' @export
conect <- function(...)
  UseMethod("conect")

#' @rdname conect
#' @export
conect.default <- function(eleid.1, eleid.2, ...)
{
  if(is.null(eleid.1) & is.null(eleid.2))
    return(NULL)
  eleid.1 <- as.integer(eleid.1)
  eleid.2 <- as.integer(eleid.2)
  con <- data.frame(eleid.1, eleid.2)
#   con[con$eleid.1 > con$eleid.2,] <- rev(con[con$eleid.1 > con$eleid.2,])
  con <- con[order(con$eleid.2),]
  con <- con[order(con$eleid.1),]
  rownames(con) <- NULL
  class(con) <- c("conect","data.frame")
  if(nrow(con) == 0)
    con <- NULL
  return(con)
}

#' @rdname conect
#' @export
conect.coords <- function(x, radii = 0.75, safety = 1.2, by.block = FALSE, ...) {
  if(!is.coords(x))
    stop("'x' must be an object of class 'coords'")

  radii <- radii*safety
  data <- cbind(x, radii)

  findCon <- function(data) {
    nat <- nrow(data)
    if(nat==0) return(NULL)
    r <- sqrt(
        outer(data$x1, data$x1, "-")^2 +
        outer(data$x2, data$x2, "-")^2 +
        outer(data$x3, data$x3, "-")^2
    )
    bond.dist <- outer(data$radii, data$radii, "+")
    M <- lower.tri(r) & (r < bond.dist)
    if(all(!M)) return(NULL)
    eleid <- matrix(rownames(data), nrow = nat, ncol = nat)   
    eleid.1 <- as.integer(t(eleid)[M])
    eleid.2 <- as.integer(  eleid [M]) 
    return(conect.default(eleid.1, eleid.2))
  }

  if(!by.block) {
    con <- findCon(data)
  }
  else {
    get.con <- function(shift = c(0,0,0), x, radii, step) {
      coords.range <- range(x)
      coords.range <- t(t(coords.range) - shift)

      x.cuts <- seq(coords.range["min","x"], coords.range["max","x"] + step, step)
      y.cuts <- seq(coords.range["min","y"], coords.range["max","y"] + step, step)
      z.cuts <- seq(coords.range["min","z"], coords.range["max","z"] + step, step)
    
      x.index <- cut(x$x1, x.cuts, include.lowest = TRUE)
      y.index <- cut(x$x2, y.cuts, include.lowest = TRUE)
      z.index <- cut(x$x3, z.cuts, include.lowest = TRUE)
      
      f <- interaction(x.index, y.index, z.index)
      
      data <- cbind(x, radii)
      
      con <- split(data, f)
      con <- lapply(con, findCon)
      con <- do.call(rbind, con)

      return(con)
    }

    width <-10
    shift <- expand.grid(0:1,0:1,0:1)*width/2

    con <- apply(shift, 1, get.con, x, radii, width)
    con <- unique(do.call(rbind, con))
    conect(con$eleid.1, con$eleid.2)
    rownames(con) <- NULL
  }
  
  return(con)
}

#' @rdname conect
#' @export
conect.pdb <- function(x, safety = 1.2, by.block = FALSE, ...) {
  symb <- toSymbols(x$atom$elename)
  symb[is.na(symb)] <- "Xx"
  rcov <- Rpdb::elements[match(symb, Rpdb::elements[,"symb"]), "rcov"]
  con <- conect(coords(x), rcov, safety, by.block)
  con <- conect(x$atoms$eleid[con$eleid.1], x$atoms$eleid[con$eleid.2])
  return(con)
}

#' @rdname conect
#' @export
is.conect <- function(x)
{
  to.return <- any(attr(x,which="class") == "conect")
  return(to.return)
}
