#' Universal Constants
#' 
#' This data set provides various universal constants
#' 
#' @format 
#' A data frame containing for each universal constant the following information.
#' \describe{
#'   \item{\code{Quantity}}{a character vector containing a short description of the constants.}
#'   \item{\code{Value}}{a numeric vector containing the value of the constants.}
#'   \item{\code{Unit}}{a character vector indicating the unit of the constants.}
#' }
#' 
#' @source 
#' http://www.ebyte.it/library/educards/constants/ConstantsOfPhysicsAndMath.html
#' 
#' @examples 
#' # Data for the speed of light
#' universalConstants["c",]
#' 
#' # Return the speed of light in m.s-1
#' universalConstants["c","Value"]
#' 
#' # Return the Planck constant in J.s
#' universalConstants["h","Value"]
#' 
#' @keywords datasets
#' 
#' @name universalConstants
NULL
