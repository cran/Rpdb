#' PDB File Writer
#' 
#' Writes a Protein Data Bank (PDB) coordinate file from an object of class \sQuote{pdb}.
#' 
#' All data stored in the \sQuote{pdb} object are written to a PDB file. A list of objects of class \sQuote{pdb} can be provided to write multiple MODEL data into a single file. In this case, each \sQuote{pdb} object of the list must have the same \code{crystal} and \code{conect} components.
#' \cr
#' To write only a subset of a \sQuote{pdb} object see function \code{\link{subset.pdb}}.
#' 
#' @param x an object, or a list of objects, of class \sQuote{pdb}.
#' @param file a single element character vector containing the name of the PDB file to be created.
#' 
#' @return No return value, called for side effects.
#'
#' @references
#' PDB format is described at:
#' http://www.wwpdb.org/documentation/format33/v3.3.html
#' 
#' @seealso \code{\link{read.pdb}}, \code{\link{pdb}}, \code{\link{crystal}}, \code{\link{atoms}}, \code{\link{conect}}, \code{\link{subset.pdb}}
#' 
#' @examples 
#' # Read a PDB file included with the package
#' pdb <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"))
#' 
#' # Write the pdb object to file "Rpdb.pdb" in the current directory
#' write.pdb(pdb, file = "Rpdb.pdb")
#'
#' # Cleanup
#' unlink("Rpdb.pdb")
#' 
#' @keywords IO
#' 
#' @name write.pdb
#' @export
write.pdb <- function(x, file="Rpdb.pdb")
{
  if(is.pdb(x)) x <- list(x)
  if( ! all(sapply(x, is.pdb)))
      stop("'x' must be an object of class 'pdb'");

  lines <- NULL;
  
  ### Title
  title = unique(unlist(lapply(x, function(y) return(y$title))));
  title = format.pdb.title.character(title);
  lines = c(lines, title);

  ### Remarks
  remark <- unique(unlist(lapply(x, function(y) return(y$remark))))
  if( ! is.null(remark))
  {
    noHeader = (substr(remark, 1, 6) != "REMARK");
    remark[noHeader] <- paste0("REMARK", remark[noHeader]);
    lines <- c(lines, remark)
  }

	### Crystal
	if(! is.null(x[[1]]$crystal))
	{
		cryst = x[[1]]$crystal;
		abc   = paste(format(cryst$abc, width=9, nsmall=3, justify="right"), collapse="")
		abg   = paste(format(cryst$abg, width=7, nsmall=2, justify="right"), collapse="")
		lines = c(lines, paste("CRYST1", abc, abg, sep=""));
	}

  for(model in seq_along(x)) {
    if(length(x) > 1) lines <- c(lines, paste0("MODEL ",format(model, width = 4)))

    if(basis(x[[model]]) != "xyz") x[[model]] <- abc2xyz.pdb(x[[model]])
  
    X <- round(x[[model]]$atoms$x1, digits = 3)
    Y <- round(x[[model]]$atoms$x2, digits = 3)
    Z <- round(x[[model]]$atoms$x3, digits = 3)
    occ  <- round(x[[model]]$atoms$occ , digits = 2)
    temp <- round(x[[model]]$atoms$temp, digits = 2)
  
    recname <- format(x[[model]]$atoms$recname, justify = "left" , width = 6)
    eleid   <- format(x[[model]]$atoms$eleid  , justify = "right", width = 5)
    elename <- format(x[[model]]$atoms$elename, justify = "left" , width = 4)
    alt     <- format(x[[model]]$atoms$alt    , justify = "left" , width = 1)
    resname <- format(x[[model]]$atoms$resname, justify = "left" , width = 3)
    chainid <- format(x[[model]]$atoms$chainid, justify = "left" , width = 1)
    resid   <- format(x[[model]]$atoms$resid  , justify = "right", width = 4)
    insert  <- format(x[[model]]$atoms$insert , justify = "left" , width = 1)
    X       <- format(X                       , justify = "right", width = 8)
    Y       <- format(Y                       , justify = "right", width = 8)
    Z       <- format(Z                       , justify = "right", width = 8)
    occ     <- format(occ                     , justify = "right", width = 6)
    temp    <- format(temp                    , justify = "right", width = 6)
    segid   <- format(x[[model]]$atoms$segid  , justify = "left" , width = 3)
  
    lines <- c(lines,paste(recname,eleid," ",elename,alt,resname," ",chainid,resid,insert,"   ",X,Y,Z,occ,temp,"      ",segid,sep=""))
    if(length(x) > 1) lines <- c(lines,"ENDMDL")
  }
  
  if(!is.null(x[[1]]$conect))
  {
    eleid.1 <- format(x[[1]]$conect$eleid.1,width=5,justify="right")
    eleid.2 <- format(x[[1]]$conect$eleid.2,width=5,justify="right")
    conect <- paste("CONECT",eleid.1,eleid.2,sep="")
    lines <- c(lines,conect)    
  }
  #
  writeLines(lines, file);
}
