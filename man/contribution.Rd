\name{contribution}
\alias{contribution}
\encoding{latin1}

\title{Contribution to propagated uncertainty for each variable}

\description{
Calculates the relative "contribution" \eqn{C_i} of each variable \eqn{x_i} to the propagated uncertainty, as defined in the \emph{Expression of the Uncertainty of Measurement in Calibration, Eqn 4.2, page 9} (see 'References'). In the implementation here, the contributions are rescaled to sum up to 1.
}

\usage{
contribution(object, plot = TRUE, ...)
}

\arguments{
  \item{object}{an object obtained from \code{\link{propagate}}.}
  \item{plot}{logical. If \code{TRUE}, a barplot with the contributions is plotted.}
  \item{...}{other parameters for \code{\link{barplot}}.}  
}

\details{
\deqn{C_i = \frac{\partial f}{\partial x_i}\sigma_i}
with \eqn{C_i} = the "contribution", which is calculated by \code{object$evalGrad * sqrt(diag(object$covMat))}.
}

\value{A named vector with relative contributions for each variable and a barplot of the values, if \code{plot = TRUE}.}

\author{
Andrej-Nikolai Spiess
}

\references{
Expression of the Uncertainty of Measurement in Calibration.\cr
European Cooperation for Accreditation (EA-4/02), 1999.
}

\examples{
EXPR1 <- expression(x^y)
x <- c(5, 0.2)
y <- c(1, 0.1)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
contribution(RES1)
}   

\keyword{algebra}

