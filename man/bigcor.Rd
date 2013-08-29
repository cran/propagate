\name{bigcor}
\alias{bigcor}
\encoding{latin1}

\title{Creating very large correlation/covariance matrices}

\description{
The storage of a value in double format needs 8 bytes. When creating large correlation matrices, the amount of RAM might not suffice, giving rise to the dreaded \emph{"cannot allocate vector of size ..."} error. For example, an input matrix with 50000 columns/100 rows will result in a correlation matrix with a size of 50000 x 50000 x 8 Byte / (1024 x 1024 x 1024) = 18.63 GByte, which is still more than most standard PCs. \code{bigcor} uses the framework of the 'ff' package to store the correlation/covariance matrix in a file. The complete matrix is created by filling a large preallocated empty matrix with sub-matrices at the corresponding positions. See 'Details'. Calculation time is ~ 20s for an input matrix of 10000 x 100 (cols x rows).
}

\usage{
bigcor(x, fun = c("cor", "cov"), size = 2000, verbose = TRUE, ...)
}

\arguments{
  \item{x}{the input matrix.}
  \item{fun}{create either a \code{\link{cor}}elation or \code{\link{cov}}ariance matrix.} 
  \item{size}{the n x n block size of the submatrices. 2000 has shown to be time-effective.}
  \item{verbose}{logical. If \code{TRUE}, information is printed in the console when running.}
  \item{...}{other parameters to be passed to \code{\link{cor}} or \code{\link{cor}}.}
}

\details{
Calculates a correlation matrix \eqn{\mathbf{C}} or covariance matrix \eqn{\mathbf{\Sigma}} using the following steps:\cr
1) An input matrix \code{x} with \eqn{N} columns is split into equal size blocks (+ a possible remainder block) of size \eqn{n}: \eqn{A_1, \ldots,  A_n;  B_1, \ldots, B_n} etc. The block size can be defined by the user, \code{size = 2000} is a good value because \code{\link{cor}} can handle this quite quickly (~ 400 ms). If the matrix has 13796 columns, the split will be 2000; 2000; 2000; 2000; 2000; 2000; 1796.\cr
2) For all combinations of blocks, the \eqn{n \times n} correlation sub-matrix is calculated, so \eqn{A/A, A/B, B/B} etc.\cr
3) The sub-matrices are transferred into a preallocated \eqn{N \times N} empty matrix at the corresponding position (where the correlations would usually reside). To ensure symmetry around the diagonal, this is done twice in the upper and lower triangle.\cr

Since the resulting matrix is in 'ff' format, one has to subset to extract regions into normal \code{\link{matrix}}-like objects. See 'Examples'.
}

\value{
The corresponding correlation/covariance matrix in 'ff' format.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\url{http://rmazing.wordpress.com/2013/02/22/bigcor-large-correlation-matrices-in-r/}
}

\examples{
## Small example to prove similarity
## to standard 'cor'. We create a matrix
## by subsetting the complete 'ff' matrix.
MAT <- matrix(rnorm(70000), ncol = 700)
COR <- bigcor(MAT, size= 2000, fun = "cor")
COR <- COR[1:nrow(COR), 1:ncol(COR)]
all.equal(COR, cor(MAT)) # => TRUE

\dontrun{
## Create large matrix.
MAT <- matrix(rnorm(137500), ncol = 13750)
COR <- bigcor(MAT, size= 2000, fun = "cor")

## Extract submatrix.
SUB <- COR[1:3000, 1:3000]
}
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

