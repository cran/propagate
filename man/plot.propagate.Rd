\name{plot.propagate}
\alias{plot.propagate}

\title{Plotting function for 'propagate' objects}

\description{
Creates two different plots from 'propagate' objects:\cr
i) a histogram of the evaluated results from the multivariate simulated data, along with a density curve and 95\% confidence intervals.\cr
ii) a boxplot of the evaluated results from the multivariate simulated data, along with first- and second-order mean/s.d. and 95\% confidence intervals.
}

\usage{
\method{plot}{propagate}(x, logx = FALSE, ...)
}

\arguments{
  \item{x}{an object returned from \code{\link{propagate}}.} 
  \item{logx}{logical. Should the data be displayed on a logarithmic abscissa?}
  \item{...}{other parameters to \code{\link{hist}} or \code{\link{boxplot}}.}  
}

\value{
A plot as described above.
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
plot(RES1)
}

\keyword{models}
\keyword{nonlinear}
