#' Basic Vectorial Operations
#' 
#' Basic vectorial operations such as scalar product and vectorial product
#' 
#' 
#' @return
#' \itemize{
#' \item{\code{dotProct} return a single element numeric vector.}
#' \item{\code{vectNorm} return a single element numeric vector.}
#' \item{\code{rotVect} return a numeric vector of length 3.}
#' \item{\code{vectProct} return a numeric vector of length 3.}
#' }
#' 
#' @param U a numeric vector of length 3.
#' @param V a numeric vector of length 3.
#' @param n an integer.
#' 
#' @seealso \code{\link[base]{matmult}}
#' 
#' @examples 
#' Vx <- c(3,0,0)
#' vectNorm(Vx)
#' Vx <- Vx / vectNorm(Vx)
#' Vy <- c(0,1,0)
#' Vz <- vectProd(Vx, Vy)
#' print(Vz)
#' 
#' @keywords manip
#' 
#' @name vectorialOperations
#' @export
dotProd <- function(U,V){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  if(!(is.vector(V) & length(V)==3 & is.numeric(V)))
    stop("'V' must be a numeric vector of length 3")
  sum(U*V)
}

#' @rdname vectorialOperations
#' @export
vectNorm <- function(U){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  sqrt(dotProd(U,U))
}

#' @rdname vectorialOperations
#' @export
rotVect <- function(U, n = 1){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  if(!(length(n)==1 & round(n)==n))
    stop("'n' must be an integer")
  if(n>0)
    U <- U[c((length(U)-n+1):length(U),1:(length(U)-n))]
  if(n<0)
    U <- U[c((abs(n)+1):length(U),1:abs(n))]
  return(U)
}

#' @rdname vectorialOperations
#' @export
vectProd <- function(U,V){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  if(!(is.vector(V) & length(V)==3 & is.numeric(V)))
    stop("'V' must be a numeric vector of length 3")
  rotVect(U,-1) * rotVect(V, 1) - rotVect(U, 1) * rotVect(V,-1)
}
  

