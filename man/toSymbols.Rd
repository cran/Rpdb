\name{toSymbols}
\alias{toSymbols}
\alias{toSymbols.integer}
\alias{toSymbols.numeric}
\alias{toSymbols.character}
\title{
Atomic Symbols Converter}
\description{Converts character strings or atomic numbers into atomic symbols.}
\usage{
toSymbols(x, \dots)

\method{toSymbols}{integer}(x, \dots)

\method{toSymbols}{numeric}(x, \dots)

\method{toSymbols}{character}(x, nletters = 3, \dots)

}
\arguments{
  \item{x}{a vector to be converted into atomic symbols.}
  \item{nletters}{an integer used to truncate the character strings before convertion.}
  \item{\dots}{further arguments passed to or from other methods.}  
}
\details{
Each elements of \code{x} are converted into atomic symbols.\cr
When \code{x} is an integer (or numeric) vector, atomic number are search into the \code{elements} data set to find associated atomic symbols.\cr
When \code{x} is a character vector, \code{toSymbols} first removes all leading and trailing white spaces and numbers.
Then translates the first character of the character strings to uppercase and all the others to lowercase.
Finally, the character strings are tested for matching with element symbols provided by the \code{elements} data set.
NA are produced for no matching.
}
\value{a character vector containing atomic symbols}

\seealso{\code{\link{elements}}}

\examples{
x <- c(1:10)
toSymbols(x)

x <- c("C  "," o","h1","1h","UU","SI0","cR")
toSymbols(x)

# 'nletters' can be used to truncate the character strings before convertion, if need
toSymbols("SIL", nletters=3) # return NA
toSymbols("SIL", nletters=2) # return "Si"
toSymbols("SIL", nletters=1) # return "S"

}
\keyword{ manip }
