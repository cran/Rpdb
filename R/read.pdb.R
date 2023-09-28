#' PDB File Reader
#' 
#' Reads a Protein Data Bank (PDB) coordinate file.
#' 
#' The \code{read.pdb} function reads the TITLE, REMARK, ATOM, HETATM, CRYST1 and CONECT records from a PDB file. Three different reading modes can be used depending on the value of \code{MODEL}: 
#' \itemize{
#'   \item When \code{MODEL} is a vector of integers, MODEL sections whose serial numbers match these integers are read.
#'   \item When \code{MODEL == NULL}, all MODEL sections are read.
#'   \item When \code{MODEL == NA}, MODEL records are ignored to read all ATOM and/or HETATM records together to return a single object.
#' }
#' 
#' @return 
#' When a single MODEL section is read, this function returns an object of class  \sQuote{pdb} (a list with a \code{class} attribute equal to \code{pdb}) with the following components:
#' \item{title}{a character vector containing the TITLE records found in the PDB file.}
#' \item{remark}{a character vector containing the REMARK records found in the PDB file.}
#' \item{cryst1}{a list of class \sQuote{cryst1} containing the first CRYST1 record found in the PDB file. All others are ignored.}
#' \item{atoms}{a data.frame of class \sQuote{atoms} containing the ATOM and HETATM records found in the PDB file.}
#' \item{conect}{a data.frame of class \sQuote{conect} containing the CONECT records found in the PDB file.}
#' When multiple MODEL sections are read, a list of object of class \sQuote{pdb} is returned.
#' 
#' @param file a single element character vector containing the name of the PDB file to be read.
#' @param ATOM a single element logical vector indicating whether ATOM records have to be read.
#' @param HETATM a single element logical vector indicating whether HETATM records have to be read.
#' @param CRYST1 a single element logical vector indicating whether CRYST1 records have to be read.
#' @param CONECT a single element logical vector indicating whether CONECT records have to be read.
#' @param TITLE a single element logical vector indicating whether TITLE records have to be read.
#' @param REMARK a single element logical vector indicating whether REMARK records have to be read.
#' @param MODEL an integer vector containing the serial number of the MODEL sections to be read. Can also be equal to NULL to read all the MODEL sections or to NA to ignore MODEL records (see details).
#' 
#' @references 
#' PDB format has been taken from:
#' http://www.wwpdb.org/documentation/format33/v3.3.html
#' 
#' @seealso 
#' \code{\link{write.pdb}}, \code{\link{pdb}}, \code{\link{cryst1}}, \code{\link{atoms}}, \code{\link{conect}}
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
read.pdb <- function(file, ATOM = TRUE, HETATM = TRUE, CRYST1 = TRUE,
                     CONECT = TRUE, TITLE = TRUE, REMARK = TRUE, MODEL = 1)
{
  if(!file.exists(file))
    stop("File '", file, "'' is missing")
  
  lines <- readLines(file)

  recname <- substr(lines, 1, 6)
  
  title <- NULL
  if(TITLE & any(recname == "TITLE "))
    title <- subset(lines, recname == "TITLE ")    

  remark <- NULL
  if(REMARK & any(recname == "REMARK"))
    remark <- subset(lines, recname == "REMARK")

  atom <- character(0)
  is.hetatm <- rep(FALSE, length(recname))
  is.atom   <- rep(FALSE, length(recname))
  if(ATOM  ) is.atom   <- recname == "ATOM  "
  if(HETATM) is.hetatm <- recname == "HETATM"
  atoms <- subset(lines, is.atom | is.hetatm)
  if(length(atoms)==0) stop("No atoms have been selected")

  model.factor <- rep(0, length(recname))
  model.ids <- "MODEL.1"
  model.start <- grep("MODEL ", recname)
  model.end   <- grep("ENDMDL", recname)
  if(length(model.start) != length(model.end))
    stop("'Unterminated MODEL section'")
  if(!(length(model.start) == 0 & length(model.end) == 0)) {
    if(any(model.start >= model.end))
      stop("'Unterminated MODEL section'")
    if(is.null(MODEL)) {
      model.ids <- as.integer(substr(lines[model.start],11,14))
      MODEL <- model.ids 
    }
    if(!is.na(MODEL[1])) {
      model.ids   <- as.integer(substr(lines[model.start],11,14))
      model.start <- model.start[model.ids %in% MODEL]
      model.end   <- model.end  [model.ids %in% MODEL]
      model.ids   <- model.ids  [model.ids %in% MODEL]
      model.ids   <- paste("MODEL",model.ids,sep=".")
      model.factor2 <- model.factor
      model.factor [model.start+1] <- model.start
      model.factor2[model.end    ] <- model.start
      model.factor <- cumsum(model.factor) - cumsum(model.factor2) ; rm(model.factor2)
      model.factor[model.factor == 0] <- NA
    }
  }
  model.factor <- as.factor(model.factor)
  levels(model.factor) <- model.ids
  model.factor <- model.factor[is.atom | is.hetatm]
  
  trim <- function(str)
    sub(' +$','',sub('^ +', '', str))
  
  recname <- trim(substr(atoms,  1,  6))
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

  atoms <- atoms.default(recname, eleid, elename, alt,
                      resname, chainid, resid, insert,
                      x1, x2, x3, occ, temp, segid, basis = "xyz")

  cryst1 <- NULL
  if(CRYST1 & any(grepl("CRYST1", lines)))
  {
    cryst1 <- subset(lines, grepl("CRYST1", lines))
    if(length(cryst1) != 1)
    {
      warning("Multiple 'CRYST1' records have been found. Only the first record has been kept.")
      cryst1 <- cryst1[1]
    }
    abc <- c(
      a = as.numeric(substr(cryst1,  7, 15)),
      b = as.numeric(substr(cryst1, 16, 24)),
      c = as.numeric(substr(cryst1, 25, 33))
    )
    abg <- c(
      alpha = as.numeric(substr(cryst1, 34, 40)),
      beta  = as.numeric(substr(cryst1, 41, 47)),
      gamma = as.numeric(substr(cryst1, 48, 54))
    )
    sgroup <- substr(cryst1, 56, 66)
    
    if(any(is.na(abc))) warning("In 'cryst1': 'abc' contains NA values")
    if(any(is.na(abg))) warning("In 'cryst1': 'abg' contains NA values")
    if(sgroup == "") sgroup <- NULL

    cryst1 <- cryst1.default(abc, abg, sgroup)
  }
  
  conect = NULL
  if(CONECT & any(grepl("CONECT", lines)))
  {
    conect <- subset(lines, grepl("CONECT", lines))
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
    conect <- subset(conect,tokeep)
  }

  to.return <- pdb(atoms, cryst1 , conect , remark , title )
  to.return <- split(to.return, model.factor)
  if(length(to.return) == 1) to.return <- to.return[[1]]

  return(to.return)
}