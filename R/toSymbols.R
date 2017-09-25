#' Atomic Symbols Converter
#' 
#' Converts character strings or atomic numbers into atomic symbols.
#' 
#' Each elements of \code{x} are converted into atomic symbols.\cr
#' When \code{x} is an integer (or numeric) vector, atomic number are search into the \code{elements} data set to find associated atomic symbols.\cr
#' When \code{x} is a character vector, \code{toSymbols} first removes all leading and trailing white spaces and numbers.
#' Then translates the first character of the character strings to uppercase and all the others to lowercase.
#' Finally, the character strings are tested for matching with element symbols provided by the \code{elements} data set.
#' NA are produced for no matching.
#' 
#' @return a character vector containing atomic symbols
#' 
#' @param x a vector to be converted into atomic symbols.
#' @param nletters an integer used to truncate the character strings before convertion.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{elements}}
#' 
#' @examples 
#' x <- c(1:10)
#' toSymbols(x)
#' 
#' x <- c("C  "," o","h1","1h","UU","SI0","cR")
#' toSymbols(x)
#' 
#' # 'nletters' can be used to truncate the character strings before convertion, if need
#' toSymbols("SIL", nletters=3) # return NA
#' toSymbols("SIL", nletters=2) # return "Si"
#' toSymbols("SIL", nletters=1) # return "S"
#' 
#' @keywords manip
#' 
#' @name toSymbols
#' @export
toSymbols <- function(x, ...)
  UseMethod("toSymbols")

#' @rdname toSymbols
#' @export
toSymbols.integer <- function(x, ...){
  if(!is.integer(x))
    stop("'x' must be an integer")
  
  symb <- Rpdb::elements[match(x, Rpdb::elements[, "num"]), "symb"]
  symb <- as.character(symb)
  
  if(any(is.na(x)))
    warning("NAs introduced by coercion")
  
  return(symb)
}

#' @rdname toSymbols
#' @export
toSymbols.numeric <- function(x, ...){
  if(!is.numeric(x))
    stop("'x' must be an numeric")
  
  if(any(round(x) != x))
    stop("'x' must be a whole number")
  
  symb <- toSymbols(as.integer(x))
  
  return(symb)
}

#' @rdname toSymbols
#' @export
toSymbols.character <- function(x, nletters = 3, ...)
{
  x <- sub(' +$','',sub('^ +', '', x))
  x <- gsub("[0-9]","",x)

  l1 <- substr(x, 1,1)
  ln <- substr(x, 2, nletters)
  
  x <- paste0(toupper(l1), tolower(ln))

  x <- Rpdb::elements[match(x, Rpdb::elements[,"symb"]),"symb"]
  x <- as.character(x)

  if(any(is.na(x)))
    warning("NAs introduced by coercion")
  
  return(x)
}