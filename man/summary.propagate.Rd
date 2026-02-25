\name{summary.propagate}
\alias{summary.propagate}

\title{Summary function for 'propagate' objects}

\description{
Provides a printed summary of the results obtained by \code{\link{propagate}}, with many GUM (2008) and GUM Supplement 1 (2008) relevant statistics, as outlined in 'Details'. 
}

\usage{
\method{summary}{propagate}(object, normality = FALSE, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{propagate}}.} 
  \item{normality}{logical. If \code{TRUE}, the three repeated subsampling normality tests are performed. See 'Details'.}
  \item{...}{other parameters for future methods.}       
}

\value{
A printed output with the items listed under 'Details'.
}

\details{
\code{summary.propagate} provides an extensive output and overview on the results of \code{\link{propagate}}:\cr
\enumerate{
\bold{For classical JCGM 100:2008 output:}
  \item Results of the first- (mean.1, u.1) and second-order (mean.2, u.2) uncertainty propagation, as well as from MC simulation (mean.MC, u.MC, median.MC, mad.MC), including confidence intervals
  \item The Welch-Satterthwaite degrees of freedom \eqn{\mathrm{ws.df}}, as obtained from \code{\link{WelchSatter}}
  \item The coverage factor \eqn{k = t_{\mathrm{ws.df}}(1-\alpha/2)}
  \item The expanded uncertainty \eqn{U = k \cdot u}, where \eqn{u} can be either first-order u.1 or second-order u.2, see \code{\link{propagate}}
  \item The symbolic gradient matrix \eqn{\nabla}, also often termed \eqn{\mathbf{G}}, \eqn{\mathbf{J}} or \eqn{\mathbf{C_x}} (JCGM 102:2011)
  \item The evaluated gradient matrix, also known as "sensitivities", \eqn{c_i = \mathrm{eval}(\partial f / \partial x_i)}  
  \item The relative contribution (summing up to 1) \deqn{\mathrm{RC}_{ii} = \mathrm{diag}\left(\frac{|c_i \cdot c_j \cdot \Sigma_{ij}|}{\sum_{k,\ell} |c_k \cdot c_\ell \cdot \Sigma_{k\ell}|}\right)}
  \item The symbolic Hessian matrix \eqn{\mathbf{H}}, also often termed \eqn{\nabla^2}
  \item The evaluated Hessian matrix \eqn{h_{ij} = \mathrm{eval}\left(\frac{\partial^2 f}{\partial x_i \partial x_j}\right)}
  \item The input covariance matrix \eqn{\mathbf{\Sigma}}, as provided from (i) the uncertainies in the 2nd row of the input \code{data}, (ii) estimation from the raw data, or (iii) an externally supplied one. Also often termed \eqn{\mathbf{V}} or \eqn{\mathbf{U_x}} (JCGM 102:2011)\cr\cr 
\bold{For JCGM 101:2008, but extended with checks of input parameters for the Copula \emph{vs.} the output attributes of the MC-generated random matrix \eqn{\mathbf{X}^\text{MC}}:}
  \item The covariance estimated from the MC simulations of the Copula distribution, \eqn{\hat{\mathbf{\Sigma}}^\text{MC} = \mathrm{cov}(\mathbf{X}^\text{MC})}
  \item The Frobenius norm of the estimated \emph{vs.} input covariance matrix, \eqn{\|\hat{\mathbf{\Sigma}}^\text{MC} - \mathbf{\Sigma} \|_F}
  \item The correlation matrix \eqn{\mathbf{R}} derived from \eqn{\mathbf{\Sigma}} (by \code{\link{cov2cor}}) to construct the Gaussian Copula with \emph{t}-margins, \eqn{\mathcal{C}^{\text{Gauss, t}}} 
  \item The correlation matrix estimated from the MC simulations of the Copula distribution, \eqn{\hat{\mathbf{R}}^\text{MC} = \mathrm{cor}(\mathbf{X}^\text{MC})}
  \item The Frobenius norm of the estimated \emph{vs.} input correlation matrix, \eqn{\|\hat{\mathbf{R}}_\text{MC} - \mathbf{R} \|_F}
  \item A comparison of Copula-input and MC-derived moments/DOFs, to check for similarity and success of the Copula construction/MC sampling\cr\cr
\bold{Further repeated subsampling tests for normality on the evaluated MC samples, \eqn{Y = \mathbf{X}^\text{MC}}:}
  \item Skewness and Excess Kurtosis 
  \item Shapiro-Wilk test (2000 repeats, 5000 samples)
  \item Anderson-Darling test (1000 repeats, 10000 samples)
  \item Lilliefors test (1000 repeats, 10000 samples)
For the three tests, the median of the log \emph{p}-values is calculated, as the resulting \emph{p}-distribution is more or less log-normal, \eqn{\mathrm{med}(P) = \mathrm{exp}(\mathrm{med}(\mathrm{log}(P)))}.    
}
}

\references{
See \code{\link{propagate}}
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
                  nsim = 10000, check = TRUE)
summary(RES1)
}

\keyword{models}
\keyword{nonlinear}
