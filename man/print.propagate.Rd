\name{print.propagate}
\alias{print.propagate}

\title{Printing function for 'propagate' objects}

\description{
Provides a printed summary of the results obtained by \code{\link{propagate}}, such as statistics of the first/second-order uncertainty propagation, Monte Carlo simulation, the covariance matrix and symbolic as well as evaluated versions of the Gradient and Hessian matrices.
}

\usage{
\method{print}{propagate}(x, verbose = TRUE,  ...)
}

\arguments{
  \item{x}{an object returned from \code{\link{propagate}}.}     
  \item{verbose}{logical. If \code{TRUE} (default), covariance, Gradient and Hessian matrices are printed.} 
  \item{...}{other parameters for future methods.}       
}

\value{
A printed summary containing the results from error propagation and Monte Carlo simulation, as well as the covariance matrix and symbolic and evaluated versions of the Gradient and Hessian matrices.
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
EXPR1 <- expression(x^2 * sin(y))
x <- c(5, 0.01)
y <- c(1, 0.01)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
print(RES1)
}

\keyword{models}
\keyword{nonlinear}
