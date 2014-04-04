\name{print.propagate}
\alias{print.propagate}

\title{Printing function for 'propagate' objects}

\description{
Prints the statistical results obtained from \code{\link{propagate}}.
}

\usage{
\method{print}{propagate}(x, ...)
}

\arguments{
  \item{x}{an object returned from \code{\link{propagate}}.}  
  \item{...}{other parameters for future methods.}       
}

\value{
A printed summary containing the results from error propagation and Monte Carlo simulation.
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
RES1
}

\keyword{models}
\keyword{nonlinear}
