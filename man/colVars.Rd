\name{matrixStats}
\alias{colVars}
\alias{rowVars}
\alias{colSDs}
\alias{rowSDs}
\encoding{latin1}

\title{fast column- and row-wise versions of variance and standard deviation for matrices}

\description{
These four functions are fast versions for column- and row-wise \code{\link{var}} and \code{\link{sd}} on matrices and are meant to substitute the classical \code{apply(mat, 1, ...)} approach.
}

\usage{
  colVars(x, na.rm = TRUE)
  rowVars(x, na.rm = TRUE)
  colSDs(x, na.rm = TRUE)
  rowSDs(x, na.rm = TRUE)  
}

\arguments{
  \item{x}{a matrix or data.frame}
  \item{na.rm}{logical. Should missing values (including NaN) be omitted from the calculations?}
}

\details{
Makes use of the derivation \deqn{Var(X)=E(X^2-2XE(X)+(E(X))^2)=E(X^2)-2(E(X))^2+(E(X))^2=E(X^2)-(E(X))^2}
}

\value{
A vector with the corresponding values.
}

\author{
  Andrej-Nikolai Spiess
}   

\examples{
## Speed comparison on large matrix.
\dontrun{
MAT <- matrix(rnorm(10 * 500000), ncol = 10)
pt1 <- proc.time()
RES1 <- apply(MAT, 1, var)
proc.time()-pt1

pt2 <- proc.time()
RES2 <- rowVars(MAT)
proc.time()-pt2

all.equal(RES1, RES2)
}
## ~ 45x speed increase
}  
   
\keyword{univar}

