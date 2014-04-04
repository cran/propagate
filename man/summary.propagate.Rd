\name{summary.propagate}
\alias{summary.propagate}

\title{Summary function for 'propagate' objects}

\description{
Provides a printed summary of the results obtained by \code{\link{propagate}}, such as statistics of the first/second-order uncertainty propagation, Monte Carlo simulation, the covariance matrix and symbolic as well as evaluated versions of the Gradient and Hessian matrices. If \code{do.sim = TRUE} had been set in \code{propagate}, skewness/kurtosis and Shapiro-Wilks/Kolmogorov-Smirnov tests for normality are calculated on the Monte-Carlo evaluations. 
}

\usage{
\method{summary}{propagate}(object, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{propagate}}.}     
  \item{...}{other parameters for future methods.}       
}

\value{
A printed summary as described above.
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
summary(RES1)
}

\keyword{models}
\keyword{nonlinear}
