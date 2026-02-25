\name{bigcor}
\alias{bigcor}
\encoding{latin1}

\title{Creating very large correlation/covariance matrices}

\description{
The storage of a value in double format needs 8 bytes. When creating large correlation matrices, the amount of RAM might not suffice, giving rise to the dreaded \emph{"cannot allocate vector of size ..."} error. For example, an input matrix with 50000 columns/100 rows will result in a correlation matrix with a size of 50000 x 50000 x 8 Byte / (1024 x 1024 x 1024) = 18.63 GByte, which still challenges most standard PCs. \code{bigcor} uses the framework of the 'hdf5r' package to store the correlation/covariance matrix in a file. The complete matrix is created by filling a large preallocated empty matrix with sub-matrices at the corresponding positions. See 'Details'.
}

\usage{
bigcor(x, y = NULL, fun = c("cor", "cov"), size = 2000, 
       verbose = TRUE, file = "result.h5", ...)
}

\arguments{
  \item{x}{the input matrix.}
  \item{y}{\code{NULL} (default) or a vector, matrix or data frame with compatible dimensions to \code{x}.}
  \item{fun}{create either a \code{\link{cor}}elation or \code{\link{cov}}ariance matrix.} 
  \item{size}{the n x n block size of the submatrices. 2000 has shown to be time-effective.}
  \item{verbose}{logical. If \code{TRUE}, information is printed in the console when running.}
  \item{file}{name (and possible path) of the created h5-file on disk.}
  \item{...}{other parameters to be passed to \code{\link{cor}} or \code{\link{cov}}.}
}

\details{
Calculates a correlation matrix \eqn{\mathbf{R}} or covariance matrix \eqn{\mathbf{\Sigma}} using the following steps:\cr
1) An input matrix \eqn{\mathbf{X}} with \eqn{N} columns is split into \eqn{k} equal size blocks (+ a possible remainder block) \eqn{A_1, A_2, \ldots,  A_k} of size \eqn{n}. The block size can be defined by the user, \code{size = 2000} is a good value because \code{\link{cor}} can handle this quite quickly (~ 400 ms). For example, if the matrix has 13796 columns, the split will be \eqn{A_1 = 1 \ldots 2000; A_2 = 2001 \ldots 4000; A_3 = 4001 \ldots 6000; A_4 = 6000 \ldots 8000 ; A_5 = 8001 \ldots 10000; A_6 = 10001 \ldots 12000; A_7 = 12001 \ldots 13796}.\cr
2) For all pairwise combinations of blocks \eqn{k \choose 2}, the \eqn{n \times n} correlation sub-matrix is calculated. If \code{y = NULL}, \eqn{\mathrm{cor}(A_1, A_1), \mathrm{cor}(A_1, A_2), \ldots, \mathrm{cor}(A_k, A_k)}, otherwise \eqn{\mathrm{cor}(A_1, y), \mathrm{cor}(A_2, y), \ldots, \mathrm{cor}(A_k, y)}.\cr
3) The sub-matrices are transferred into a preallocated \eqn{N \times N} empty matrix at the corresponding position (where the correlations would usually reside). To ensure symmetry around the diagonal, this is done twice in the upper and lower triangle. If \code{y} was supplied, a \eqn{N \times M} matrix is filled, with \eqn{M} = number of columns in \code{y}.

Since the resulting matrix is in 'hdf5r' format, one has to subset to extract regions into normal \code{\link{matrix}}-like objects. See 'Examples'.
}

\value{
The corresponding correlation/covariance matrix in 'hdf5r' format.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\url{https://rmazing.wordpress.com/2013/02/22/bigcor-large-correlation-matrices-in-r/}
}

\examples{
\donttest{
## Small example to prove similarity
## to standard 'cor'. We create a matrix
## by subsetting the complete 'hdf5r' matrix.
MAT <- matrix(rnorm(70000), ncol = 700)
h5file <- tempfile(fileext = ".h5")
COR1 <- bigcor(MAT, size = 500, fun = "cor", file = h5file)
h5 <- hdf5r::H5File$new(COR1, mode = "r")
COR1 <- h5[["matrix"]][,]  # Read entire matrix
## COR1 <- h5[["matrix"]][1:1000, 1:1000]  # would be submatrix
h5$close_all()
all.equal(COR1, cor(MAT)) # => TRUE
file.remove(h5file)
}
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

