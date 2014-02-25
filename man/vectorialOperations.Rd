\name{vectorialOperations}
\alias{dotProd}
\alias{vectNorm}
\alias{rotVect}
\alias{vectProd}
\title{Basic Vectorial Operations}
\description{Basic vectorial operations such as scalar product and vectorial product}
\usage{
## Scalar product
dotProd(U, V)

## Norm of a vector
vectNorm(U)

## Rotation of the components of a vector
rotVect(U, n = 1)

## Vectorial product
vectProd(U, V)
}
\arguments{
  \item{U}{a numeric vector of length 3.}
  \item{V}{a numeric vector of length 3.}
  \item{n}{an integer.}
}

\value{
\itemize{
\item{\code{dotProct} return a single element numeric vector.}
\item{\code{vectNorm} return a single element numeric vector.}
\item{\code{rotVect} return a numeric vector of length 3.}
\item{\code{vectProct} return a numeric vector of length 3.}
}
}

\seealso{
\code{\link{matmult}}
}
\examples{
Vx <- c(3,0,0)
vectNorm(Vx)
Vx <- Vx/vectNorm(Vx)
Vy <- c(0,1,0)
Vz <- vectProd(Vx, Vy)
print(Vz)
}
\keyword{ manip }
