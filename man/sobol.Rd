\name{sobol}
\alias{sobol}

\title{Compute Sobol Sensitivity Indices from a \code{propagate} object}

\description{
  Compute first-order and total-order Sobol sensitivity indices from an object returned by \code{\link{propagate}}, by performing evaluations based on the MC simulated samples of the provided \code{propagate} object, or from external matrices/expressions.
}

\usage{
  sobol(prop = NULL, A = NULL, B = NULL, expr = NULL, 
        method = c("jansen", "sobol", "saltelli", "homma"))
}

\arguments{
  \item{prop}{A \code{propagate} object.}
  \item{A}{An MC simulated matrix.}
  \item{B}{A second MC simulated matrix, derived from the same distributions.}
  \item{expr}{an \code{expression} to be evaluated.}
  \item{method}{Method for sensitivity indices. See 'Details'.}
}

\value{
A list with the following components:
  \item{S}{Vector of first-order indices \eqn{S_i}.}
  \item{ST}{Vector of total-order indices \eqn{S_{T_i}}.}
}

\details{
\bold{Sobol (1993) estimator:}
\deqn{
S_i = \frac{1}{V} , \mathbb{E}\left[ f(\mathbf{B}) \left( f(\mathbf{A}_B^{(i)}) - f(\mathbf{A}) \right) \right]
}
Estimates the first-order (main effect) Sobol index, \emph{i.e.} the fraction of output variance explained by input \eqn{X_i} alone, excluding all interaction effects. This is the original Monte Carlo estimator associated with the Sobol variance decomposition.

\bold{Saltelli (2002) estimator:}
\deqn{
S_i = \frac{1}{V} , \mathbb{E}\left[ f(\mathbf{B}) \left( f(\mathbf{A}_B^{(i)}) - f(\mathbf{A}) \right) \right]
}
\deqn{
S_{T_i} = \frac{1}{2V} , \mathbb{E}\left[ \left( f(\mathbf{A}) - f(\mathbf{A}_B^{(i)}) \right)^2 \right]
}
Provides both first-order and total-order Sobol indices using a single sampling design. This formulation improves numerical efficiency over the original Sobol estimator and is widely used in practice.

\bold{Homma-Saltelli (1996) estimator:}
\deqn{
S_{T_i} = \frac{1}{V} , \mathbb{E}\left[ f(\mathbf{A}) \left( f(\mathbf{A}) - f(\mathbf{A}_B^{(i)}) \right) \right]
}
Estimates the total-order Sobol index only, \emph{i.e.} the total contribution of input \eqn{X_i} to output variance including all interaction effects. No first-order estimator is defined in this method.

\bold{Jansen (1999) estimator (recommended):}
\deqn{
S_i = 1 - \frac{1}{2V} , \mathbb{E}\left[ \left( f(\mathbf{B}) - f(\mathbf{A}_B^{(i)}) \right)^2 \right]
}
\deqn{
S_{T_i} = \frac{1}{2V} , \mathbb{E}\left[ \left( f(\mathbf{A}) - f(\mathbf{A}_B^{(i)}) \right)^2 \right]
}
Provides low-variance, numerically stable estimators for both first-order and total-order Sobol indices. This formulation avoids covariance terms and is considered best practice in modern global sensitivity analysis.
}

\references{
Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index.\cr
Saltelli A, Annoni P, Azzini I, Campolongo F, Ratto M, Tarantola S.\cr
\emph{Comp Phys Comm} (2010), \bold{181}: 259-270.

Making best use of model evaluations to compute sensitivity indices.\cr
Saltelli A.\cr
\emph{Comp Phys Comm} (2002), \bold{145}: 280-297.

Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates.\cr
Sobol IM.\cr
\emph{Math Comp Sim} (2001), \bold{55}: 271-280.

Importance measures in global sensitivity analysis of nonlinear models.\cr
Homma T, Saltelli A.\cr
\emph{Rel Engineer & Sys Safety} (1996), \bold{52}: 1-17.

Analysis of variance designs for model output.\cr
Jansen MJW.\cr
\emph{Comp Phys Comm} (1999), \bold{117}: 35-43.

An importance quantification technique in uncertainty analysis for computer models.\cr
Ishigami T, Homma T.\cr
Japan Atomic Energy Research Inst Tokyo (1989).
}

\examples{
## Example from 'propagate'
EXPR1 <- expression(x/y)
x <- c(5, 0.01, 12)
y <- c(1, 0.01, 5)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", nsim = 10000, check = FALSE)
sobol(RES1)

## Classical Ishigami function for gauging
A <- cbind(
  x1 = runif(1e6, -pi, pi),
  x2 = runif(1e6, -pi, pi),
  x3 = runif(1e6, -pi, pi)
)
B <- cbind(
  x1 = runif(1e6, -pi, pi),
  x2 = runif(1e6, -pi, pi),
  x3 = runif(1e6, -pi, pi)
)
EXPR <- expression(sin(x1) + 7 * (sin(x2))^2 + 0.1 * x3^4 * sin(x1))
sobol(A = A, B = B, expr = EXPR, method = "jansen")
# => X1: 0.559, X2: 0.442, X3: 0.242
# => Homma & Saltelli (1996), Table 4: X1: 0.557, X2: 0.444, X3: 0.241
}

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

