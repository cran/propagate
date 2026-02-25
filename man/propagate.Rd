\name{propagate}
\alias{propagate}
\encoding{latin1}

\title{Propagation of uncertainty using higher-order Taylor expansion and Monte Carlo simulation}

\description{
A general function for the calculation of uncertainty propagation by first-/second-order Taylor expansion and Monte Carlo (MC) simulation including covariances. Input data can be any symbolic/numeric differentiable expression and data based on summaries (Mean/SEM/DOF), experimental measurements or sampled from distributions. Uncertainty propagation is based completely on matrix calculus accounting for full covariance structure. Monte Carlo simulation is conducted using a Gaussian copula with marginal \emph{t}-distributions. Propagation confidence intervals are calculated from the expanded uncertainties by \code{\link{WelchSatter}}'s total degrees-of-freedom or from the \eqn{\alpha}-quantiles of the MC evaluations.
}

\usage{
propagate(expr, data, type = c("stat", "raw", "sim"), cov = c("diag", "full"), 
          df.tot = NULL, k = NULL, nsim = 1000000, alpha = 0.05, 
          exp.uc = c("1st", "2nd"), check = TRUE, ...)  
}

\arguments{
  \item{expr}{an expression, such as \code{expression(x/y)}.}
  \item{data}{a dataframe or matrix containing either a) the means (\eqn{1^{st}} row), standard error of the means (\eqn{2^{nd}} row) and optional degrees of freedom (\eqn{3^{rd}} row), b) raw measurement values or c) sampled distribution data generated from any of \R's \code{\link{distributions}} or those implemented in this package (\code{\link{rDistr}}).} 
  \item{type}{type of input data. Either \code{stat} for summary data, \code{raw} for experimental measurements or \code{sim} for distribution data. See 'Details'.}  
  \item{cov}{use only \code{"diag"}onal variances, the \code{"full"} covariance matrix with off-diagonals or an external covariance matrix with the same column names as \code{data}. See 'Details'.}
  \item{df.tot}{an optional total degrees of freedom \eqn{\nu_{\mathrm{tot}}} that overrides the one obtained from \code{\link{WelchSatter}}.}
  \item{k}{an optional coverage factor that overrides the one obtained from \code{\link{WelchSatter}}.}
  \item{nsim}{the number of Monte Carlo simulations to be performed, minimum is 10000.}  
  \item{alpha}{the 1 - confidence level.}
  \item{exp.uc}{Should the expanded uncertainty and confidence interval be calculated from the \code{"1st"}-order or \code{"2nd"}-order uncertainty?}
  \item{check}{logical. If \code{TRUE}, the Copula margins are fitted with a scaled/shifted \emph{t}-distribution to check for proper construction, see \code{\link{summary.propagate}}.}
  \item{...}{other parameters to be supplied to future methods.}
}

\details{
The different \code{type} definitions are the following:\cr\cr
1) \bold{"stat", GUM 4.2: Statistical summary data (Mean, SEM, DOF), Type A evaluation}
\preformatted{
     V1            V2
|------------|------------|
|   Mean_1   |   Mean_2   |  => Means in 1st row
|------------|------------|
|   SEM_1    |   SEM_2    |  => SEMs (uncertainties) in 2nd row
|------------|------------|
|  (DOF_1)   |  (DOF_2)   |  => optional DOFs in 3rd row 
|------------|------------|
}
Covariance matrix: from \eqn{\text{SEMs}^2} as diagonals (\code{cov = "diag"}) or supplied externally.\cr

2) \bold{"raw", GUM 4.2: Raw measurement values in rows, Type A evaluation}
\preformatted{
      V1           V2
|------------|------------| -- 
|    x_1     |    y_1     |   \
|------------|------------|    \
|    x_2     |    y_2     |     |
|------------|------------|     |- Single raw values in rows 1...n 
|    x_3     |    y_3     |     |
|------------|------------|    /
|    x_n     |    y_n     |   /
|------------|------------| --
}
Mean, SEM and DOF are estimated from the data.\cr
Covariance matrix: from estimated \eqn{\text{SEMs}^2} as diagonals (\code{cov = "diag"}), full covariance matrix \eqn{\mathrm{cov}(\mathrm{data})/N} (\code{cov = "full"}) or supplied externally.\cr

3) \bold{"sim", GUM 4.3: Random samples from distributions \eqn{x_i, y_i, ... \sim F_1, F_2, ...}, Type B evaluation}
\preformatted{
      F1           F2
|------------|------------| -- 
|     x_1    |    y_1     |   \
|------------|------------|    \
|     x_2    |    y_2     |     |
|------------|------------|     |- MC samples in rows 1...n 
|     x_3    |    y_3     |     |
|------------|------------|    /
|     x_n    |    y_n     |   /
|------------|------------| --
}
Mean, SD and DOF are estimated from the data. The latter will be high (\emph{N} - 1) and approximate a normal distribution.\cr
Covariance matrix: from estimated \eqn{\text{SDs}^2} as diagonals (\code{cov = "diag"}), full covariance matrix \eqn{\mathrm{cov}(\mathrm{data})} (\code{cov = "full"}) or supplied externally.\cr

For all \code{type}s: Expanded uncertainty U and confidence intervals from \eqn{U = k \cdot u, k = t_{\mathrm{ws.df}}(1-\alpha/2), \mathrm{CI} = [\hat{y} - U, \hat{y} + U]}, and empirical \eqn{[\alpha/2, 1-\alpha/2]} quantiles of the MC output distribution.\cr\cr
\code{propagate} tries "hard" to guess the data \code{type} and gives a comment (no change in settings!), \emph{i.e.} \bold{Type A summary} => 2/3 rows, \bold{Type A raw measurement} => 3 < rows < 1000, \bold{Type B distributions} => 1000 < rows.   

The implemented methods are:\cr
1) \bold{Uncertainty propagation:}\cr
The propagated uncertainties are calculated by first-/second-order Taylor expansion accounting for full covariance structure using matrix algebra.\cr
The following transformations based on two variables \eqn{x_1, x_2} illustrate the equivalence of the matrix-based approach with well-known classical notations:\cr
\bold{First-order mean:} \eqn{\rm{E[y]} = f(\bar{x}_i)}\cr
\bold{First-order variance:} \eqn{\sigma_y^2 = {\color{red} \nabla \mathbf{\Sigma} \nabla^T}}:\cr
\deqn{{ \color{red}[\rm{j_1}\; \rm{j_2}] \left[ \begin{array}{cc} \sigma_1^2 & \sigma_1\sigma_2 \\ \sigma_2\sigma_1 & \sigma_2^2 \end{array} \right] \left[ \begin{array}{c} \rm{j_1} \\ \rm{j_2} \end{array} \right]} = \rm{j_1}^2 \sigma_1^2 + \rm{2 j_1 j_2} \sigma_1 \sigma_2 + \rm{j_2}^2 \sigma_2^2}
\deqn{= \underbrace{\sum_{i=1}^2 \rm{j_i}^2 \sigma_i^2 + 2\sum_{i=1\atop i \neq k}^2\sum_{k=1\atop k \neq i}^2 \rm{j_i j_k} \sigma_{ik}}_{\rm{classical\;notation}} = \frac{1}{1!} \left(\sum_{i=1}^2 \frac{\partial f}{\partial x_i} \sigma_i \right)^2}
\cr
\bold{Second-order mean:} \eqn{\rm{E}[y] = f(\bar{x}_i) + {\color{blue} \frac{1}{2}\rm{tr}(\mathbf{H\Sigma)}}}:\cr
\deqn{{ \color{blue} \frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \sigma_1^2 & \sigma_1\sigma_2 \\ \sigma_2\sigma_1 & \sigma_2^2 \end{array} \right]} = \frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} \sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 & \rm{h_1}\sigma_1\sigma_2 + \rm{h_2}\sigma_2^2 \\ \rm{h_3} \sigma_1^2 + \rm{h_4} \sigma_1\sigma_2 & \rm{h_3} \sigma_1\sigma_2 + \rm{h_4} \sigma_2^2 \end{array} \right]}
\deqn{ = \frac{1}{2}(\rm{h_1}\sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 + \rm{h_3}\sigma_1\sigma_2 + \rm{h_4}\sigma_2^2) = \frac{1}{2!} \left(\sum_{i=1}^2 \frac{\partial}{\partial x_i} \sigma_i \right)^2 \it f} 
\cr
\bold{Second-order variance:} \eqn{\sigma_y^2 = {\color{red} \nabla\mathbf{\Sigma}\nabla^T} + {\color{blue} \frac{1}{2}\rm{tr}(\mathbf{H\Sigma H\Sigma)}}}:\cr
\deqn{{\color{blue}\frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \rm{\sigma_1^2} & \rm{\sigma_1\sigma_2} \\ \rm{\sigma_2\sigma_1} & \rm{\sigma_2^2} \end{array} \right] \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \rm{\sigma_1^2} & \rm{\sigma_1\sigma_2} \\ \rm{\sigma_2\sigma_1} & \rm{\sigma_2^2} \end{array} \right]} = \ldots}
\deqn{= \frac{1}{2} (\rm{h_1}^2\sigma_1^4 + \rm{2h_1h_2}\sigma_1^3\sigma_2 + \rm{2h_1h_3}\sigma_1^3\sigma_2 + \rm{h_2}^2\sigma_1^2\sigma_2^2 + \rm{2h_2h_3}\sigma_1^2\sigma_2^2 + \rm{h_3}^2\sigma_1^2\sigma_2^2 + \rm{2h_1h_4}\sigma_1^2\sigma_2^2}
\deqn{+ \rm{2h_2h_4}\sigma_1\sigma_2^3 + \rm{2h_3h_4}\sigma_1\sigma_2^3 + \rm{h_4}^2\sigma_2^4 = \frac{1}{2} (\rm{h_1}\sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 + \rm{h_3}\sigma_1\sigma_2 + \rm{h_4}\sigma_2^2)^2}
\deqn{= \frac{1}{2!} \left( \left(\sum_{i=1}^2 \frac{\partial}{\partial x_i} \sigma_i \right)^2 \it f \right)^2}
\cr
with \eqn{\mathrm{E}(y)} = expectation of \eqn{y}, \eqn{\mathbf{\sigma_y^2}} = variance of \eqn{y}, \eqn{{\color{red} \nabla}} = the p x n gradient matrix with all partial first derivatives \eqn{{\color{red} \rm{j_i}}}, \eqn{\mathbf{\Sigma}} = the p x p covariance matrix, \eqn{{\color{blue}\mathbf{H}}} the Hessian matrix with all partial second derivatives \eqn{{\color{blue} \rm{h_i}}}, \eqn{\sigma_i} = the uncertainties and \eqn{\rm{tr}(\cdot)} = the trace (sum of diagonal) of a matrix. Note that because the Hessian matrices are symmetric, \eqn{{\color{blue} \rm{h_2}} = {\color{blue} \rm{h_3}}}. For a detailed derivation, see 'References'.\cr
The second-order Taylor expansion corrects for bias in nonlinear expressions as the first-order Taylor expansion assumes linearity around \eqn{\bar{X}_i}.\cr
For setups in which there is no symbolic derivation possible (i.e. \code{e <- expression(abs(x))} => "Function 'abs' is not in the derivatives table") the function automatically switches from symbolic (using \code{\link{makeGrad}} or \code{\link{makeHess}}) to numeric (\code{\link{numGrad}} or \code{\link{numHess}}) differentiation.\cr
The function will try to evaluate the expression in an environment using \code{\link{eval}} which results in a significant speed enhancement (~ 10-fold). If that fails, evaluation is done over the rows of the simulated data using \code{\link{apply}}.\cr

2) \bold{Monte Carlo simulation with a Gaussian Copula of \emph{t}-Distributed Marginals:}\cr
According to Sklar's theorem, the joint distribution of correlated random variables \eqn{X_1, \ldots, X_d} with marginal CDFs \eqn{F_1, \ldots, F_d} can be expressed as:
\deqn{F(x_1, \ldots, x_d) = C_{\mathbf{R}}^{\text{Gauss}}(F_1(x_1), \ldots, F_d(x_d))}
where \eqn{C_{\mathbf{R}}^{\text{Gauss}}} is the Gaussian copula with correlation matrix \eqn{\mathbf{R}}. For \emph{t}-distributed marginals with individual degrees of freedom \eqn{\nu_1, \ldots, \nu_d}, the marginal CDFs are:
\deqn{F_i(x_i) = F_{t,\nu_i}\left(\frac{x_i - \mu_i}{\sigma_i}\right)}
where \eqn{\mu_i} and \eqn{\sigma_i} are the location and scale parameters, and \eqn{F_{t,\nu_i}} is the CDF of the standard \emph{t}-distribution with \eqn{\nu_i} degrees of freedom. 
Sampling is as follows:
\enumerate{
  \item Convert covariance matrix \eqn{\mathbf{\Sigma}} to correlation matrix \eqn{\mathbf{R}}: \eqn{\mathbf{R}_{ij} = \mathbf{\Sigma}_{ij}/\sqrt{\mathbf{\Sigma}_{ii}\mathbf{\Sigma}_{jj}}}
  \item Generate \eqn{\mathbf{v} \sim C_{\mathbf{R}}^{\text{Gauss}}} from the Gaussian copula (uniform marginals on [0,1])
  \item Transform: \eqn{z_i = F_{t,\nu_i}^{-1}(\mathbf{v}_i)} for each margin
  \item Apply variance correction: \eqn{z_i^* = z_i / \sqrt{\nu_i/(\nu_i-2)}} for \eqn{\nu_i > 2}
  \item Scale and shift: \eqn{X_i = \mu_i + \sigma_i \cdot z_i^*} giving final MC-simulated matrix \eqn{\mathbf{X}^\text{MC}}.
}

The variance correction in step 3 is necessary because the standard \emph{t}-distribution has variance \eqn{\nu_i/(\nu_i-2)} rather than unity. 
The Gaussian copula preserves the correlation structure \eqn{\mathbf{R}} while allowing each marginal to have its own degrees of freedom \eqn{\nu_i}, making it suitable for GUM Supplement 1 uncertainty propagation with heterogeneous Type A uncertainties.\cr

The function conducts the following three additional checks, visible in \code{\link{summary.propagate}}:
\enumerate{
  \item Frobenius Norm of input covariance matrix \eqn{\mathbf{\Sigma}} to MC-derived covariance matrix \eqn{\mathbf{\hat{\Sigma}}^{\text{MC}}}: \deqn{\varepsilon_F = \frac{\sqrt{\sum_{i=1}^{n}\sum_{j=1}^{n}\left(\mathbf{\hat{\Sigma}}_{ij}^{\text{MC}}-\mathbf{\Sigma}_{ij}\right)^2}}{\sqrt{\sum_{i=1}^{n}\sum_{j=1}^{n}\mathbf{\Sigma}_{ij}^2}}}
  \item Frobenius Norm of input correlation matrix \eqn{\mathbf{R}} to MC-derived correlation matrix \eqn{\mathbf{\hat{R}}^{\text{MC}}}: \deqn{\varepsilon_F = \frac{\sqrt{\sum_{i=1}^{n}\sum_{j=1}^{n}\left(\hat{\mathbf{R}}_{ij}^{\text{MC}}-\mathbf{R}_{ij}\right)^2}}{\sqrt{\sum_{i=1}^{n}\sum_{j=1}^{n}\mathbf{R}_{ij}^2}}}
  \item For each marginal component \eqn{X_i} of the constructed copula, a shifted and scaled \emph{t}-distribution is fitted independently by estimating \eqn{\hat{\mu}_i}, \eqn{\hat{\sigma}_i}, and \eqn{\hat{\nu}_i}, checking for similarity of the margins to the input structure.
}
}

\value{
A short printed output of the results (for a more detailed output use \code{summary(object)}) and a list with the following components:   
\item{gradient}{the symbolic gradient vector \eqn{\nabla} of partial first-order derivatives.}
\item{evalGrad}{the evaluated gradient vector \eqn{\nabla} of partial first-order derivatives, also known as "sensitivities". See \code{\link{summary.propagate}}.}
\item{hessian}{the symbolic Hessian matrix \eqn{\mathbf{H}} of partial second-order derivatives.}
\item{evalHess}{the evaluated Hessian matrix \eqn{\mathbf{H}} of partial second-order derivatives.}
\item{rel.contr}{the relative contributions, \eqn{\text{RC}_i}, see \code{\link{summary.propagate}}.}
\item{covMat}{the covariance matrix \eqn{\mathbf{\Sigma}} used for Monte Carlo simulation and uncertainty propagation (input).}
\item{covMat.MC}{the estimated covariance matrix \eqn{\mathbf{\hat{\Sigma}}^{\text{MC}}} derived from the Monte Carlo simulation (output).}
\item{corMat}{the \eqn{\mathbf{\Sigma}}-converted correlation matrix \eqn{\mathbf{R}} used for Monte Carlo simulation (input).}
\item{corMat.MC}{the estimated correlation matrix \eqn{\mathbf{\hat{R}}^{\text{MC}}} derived from the Monte Carlo simulation (output).}
\item{frobCOV}{the Frobenius norm \eqn{\varepsilon_F} of input and Copula-derived covariance matrices.}
\item{frobCOR}{the Frobenius norm \eqn{\varepsilon_F} of input and Copula-derived correlation matrices.}
\item{ws.df}{the Welch-Satterthwaite degrees of freedom \eqn{\nu_{\mathrm{ws}}}, as obtained from \code{\link{WelchSatter}}.}
\item{k}{the coverage factor \eqn{k}, as calculated by \eqn{t(1-(\alpha/2), \nu_{\mathrm{ws}})}.}
\item{u.exp}{the expanded uncertainty, \eqn{U = k \cdot u}, where \eqn{u} is derived from the first- or second-order propagation, depending on the setting of \code{exp.uc}.}
\item{resSIM}{a vector containing the \code{nsim} values obtained from the expression evaluations \eqn{f(x_{m, i})} of the simulated data in \code{datSIM}.}
\item{datSIM}{a vector containing the \code{nsim} simulated multivariate values for each variable in column format.}  
\item{datSIM2}{a second Copula function MC simulation to be used as matrix \eqn{\mathbf{B}} for Sobol sensitivity analysis, see \code{\link{sobol}}.}  
\item{prop}{a summary vector containing first-/second-order expectations and uncertainties from propagation as well as the confidence interval based on the \emph{t}-distribution.} 
\item{sim}{a summary vector containing the mean, standard deviation, median, MAD from MC evaluation, as well as the confidence interval based on the \eqn{\alpha}-quantiles.}
\item{expr}{the original expression \code{expr}.}
\item{data}{a matrix containing the given (\code{type = "stat"}) or estimated (\code{type = "raw/sim"}) means, uncertainties and degrees-of-freedom.}
\item{alpha}{the original \eqn{\alpha}-level.}
\item{nsim}{the number of Monte Carlo samples used for the Copula.}
\item{check}{a matrix containing the means, uncertainties and degrees-of-freedom derived from fitting the Copula margins with a scaled/shifted \emph{t}-distribution.}
\item{type}{the original \code{type} definition, important for, \emph{e.g.}, \code{\link{sobol}}.}
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\bold{Uncertainty propagation (in general):}\cr
An Introduction to error analysis.\cr
Taylor JR.\cr
University Science Books (1996), New York.

Evaluation of measurement data - Guide to the expression of uncertainty in measurement.\cr
JCGM 100:2008 (GUM 1995 with minor corrections).\cr
\url{https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf/}.

Evaluation of measurement data - Supplement 1 to the Guide to the expression of uncertainty in measurement - Propagation of distributions using a Monte Carlo Method.\cr
JCGM 101:2008.\cr
\url{https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf/}.

\bold{Higher-order Taylor expansion:}\cr
On higher-order corrections for propagating uncertainties.\cr
Wang CM & Iyer HK.\cr
\emph{Metrologia} (2005), \bold{42}: 406-410.

Propagation of uncertainty: Expressions of second and third order uncertainty with third and fourth moments.\cr
Mekid S & Vaja D.\cr
\emph{Measurement} (2008), \bold{41}: 600-609.

\bold{Matrix algebra for uncertainty propagation:}\cr
An Introduction to Error Propagation: Derivation, Meaning and Examples of Equation Cy = FxCxFx^t.\cr
\url{https://www.research-collection.ethz.ch/handle/20.500.11850/82620}.

Second order nonlinear uncertainty modeling in strapdown integration using MEMS IMUs.\cr
Zhang M, Hol JD, Slot L, Luinge H.\cr
2011 Proceedings of the 14th International Conference on Information Fusion (FUSION) (2011).

Uncertainty propagation in non-linear measurement equations.\cr
Mana G & Pennecchi F.\cr
\emph{Metrologia} (2007), \bold{44}: 246-251.

A compact tensor algebra expression of the law of propagation of uncertainty.\cr
Bouchot C, Quilantan JLC, Ochoa JCS.\cr
\emph{Metrologia} (2011), \bold{48}: L22-L28.

Nonlinear error propagation law.\cr
Kubacek L.\cr
\emph{Appl Math} (1996), \bold{41}: 329-345.

\bold{Monte Carlo simulation (normal- and t-distribution):}\cr
MUSE: computational aspects of a GUM supplement 1 implementation.\cr
Mueller M, Wolf M, Roesslein M.\cr
\emph{Metrologia} (2008), \bold{45}: 586-594.

Copulas for uncertainty analysis.\cr
Possolo A.\cr
\emph{Metrologia} (2010), \bold{47}: 262-271.

Random variables, joint distribution functions, and copulas.\cr
Sklar A.\cr
\emph{Kybernetika} (1973), \bold{9}: 449-460.

\bold{Multivariate normal distribution:}\cr
Stochastic Simulation.\cr
Ripley BD.\cr
Stochastic Simulation (1987). Wiley. Page 98.

\bold{Testing for normal distribution:}\cr
Testing for  Normality.\cr
Thode Jr. HC.\cr
Marcel Dekker (2002), New York.

Approximating the Shapiro-Wilk W-test for non-normality.\cr
Royston P.\cr 
\emph{Stat Comp} (1992), \bold{2}: 117-119.
}

\examples{
## ATTENTION: all examples with nsim = 10000
## and check = FALSE for reduced CRAN check time.
## It is advised to use nsim = 1000000 (default)
## and check = TRUE in real-life setup.

## Example without given degrees-of-freedom.
EXPR1 <- expression(x/y)
x <- c(5, 0.01)
y <- c(1, 0.01)
DF1 <- cbind(x, y)
propagate(expr = EXPR1, data = DF1, type = "stat", 
          nsim = 10000, check = FALSE)

## Same example with given degrees-of-freedom
## => third row in input 'data'.
EXPR2 <- expression(x/y)
x <- c(5, 0.01, 12)
y <- c(1, 0.01, 5)
DF2 <- cbind(x, y)
propagate(expr = EXPR2, data = DF2, type = "stat", 
          nsim = 10000, check = FALSE)

\donttest{
################# GUM 2008 (1) ########################
## Example in Annex H.1 from the GUM 2008 manual
## (see 'References'), an end gauge calibration
## study. We use total df = 16 and alpha = 0.01, 
## as detailed in GUM H.1.6.
EXPR3 <- expression(ls + d - ls * (da * the + as * dt))
ls <- c(50000623, 25)
d <- c(215, 9.7)
da <- c(0, 0.58E-6)
the <- c(-0.1, 0.41)
as <- c(11.5E-6, 1.2E-6)
dt <- c(0, 0.029)
DF3 <- cbind(ls, d, da, the, as, dt)
propagate(expr = EXPR3, data = DF3, type = "stat", 
          nsim = 10000, check = FALSE,df.tot = 16, alpha = 0.01)
## 1st order => propagate: u.1 = 31.71, GUM H.1.4/H.6c: u = 32  
## Expanded uncertainty => propagate: 92.62, GUM H.1.6: 93
## 2nd order => propagate: u.2 = 33.91115, GUM H.1.7: u = 34.
## Also similar to the non-matrix-based approach
## in Wang et al. (2005, page 408): u1 = 33.91115.
## NOTE: After second-order correction ("u.2"), 
## uncertainty is more similar to the uncertainty
## obtained from Monte Carlo simulation!

#################### GUM 2008 (2) #################
## Example in Annex H.2 from the GUM 2008 manual
## (see 'References'), simultaneous resistance
## and reactance measurement.
data(H.2)

## This gives exactly the means, uncertainties and
## correlations as given in Table H.2:
colMeans(H.2)
sqrt(colVarsC(H.2))/sqrt(5)
cor(H.2)

## H.2.3 Approach 2 using the raw
## measurement values and full covariance
EXPR4 <- expression((V/I) *  cos(phi)) ## R
propagate(expr = EXPR4, data = H.2, type = "raw", cov = "full",
          nsim = 10000, check = FALSE)
## This gives exactly the means and uncertainties
## as given in H.2.4, Table H.4

## H.2.3 Approach 2 using the raw
## measurement values and only diagonals (no correlation)
propagate(expr = EXPR4, data = H.2, type = "raw",
          nsim = 10000, check = FALSE)
## This gives exactly the means and uncertainties
## as given in H.2.4, Table H.5

######### GUM 2008 Supplement 1 (1) #######################
## Example from 9.2.2 of the GUM 2008 Supplement 1
## (see 'References'), normally distributed input
## quantities. Assign values as in 9.2.2.1.
EXPR5 <- expression(X1 + X2 + X3 + X4)
X1 <- c(0, 1)
X2 <- c(0, 1)
X3 <- c(0, 1)
X4 <- c(0, 1)
DF5 <- cbind(X1, X2, X3, X4)
propagate(expr = EXPR5, data = DF5, type = "stat", 
          nsim = 10000, check = FALSE)
## This will give exactly the same results as in 
## 9.2.2.6, Table 2.

######### GUM 2008 Supplement 1 (2) #######################
## Example from 9.3 of the GUM 2008 Supplement 1
## (see 'References'), mass calibration.
## Formula 24 in 9.3.1.3 and values as in 9.3.1.4, Table 5.
EXPR6 <- expression((Mrc + dMrc) * (1 + (Pa - Pa0) * ((1/Pw) - (1/Pr))) - Mnom)
Mrc <- rnorm(1E5, 100000, 0.050)
dMrc <- rnorm(1E5, 1.234, 0.020)
Pa <- runif(1E5, 1.10, 1.30)  ## E(Pa) = 1.2, (b-a)/2 = 0.1 
Pw <- runif(1E5, 7000, 9000)  ## E(Pw) = 8000, (b-a)/2 = 1000
Pr <- runif(1E5, 7950, 8050) ## E(Pr) = 8000, (b-a)/2 = 50
Pa0 <- 1.2 
Mnom <- 100000
DF6 <- cbind(Mrc, dMrc, Pa, Pw, Pr, Pa0, Mnom)
propagate(expr = EXPR6, data = DF6, type = "sim", check = FALSE)
## This will give exactly the same results as in 
## 9.3.2.3, Table 6

######### GUM 2008 Supplement 1 (3) #######################
## Example from 9.4 of the GUM 2008 Supplement 1
## (see 'References'), comparison loss in microwave
## power meter calibration, zero covariance.
## Formula 28 in 9.4.1.5 and values as in 9.4.1.7.
EXPR7 <- expression(X1^2 - X2^2)
X1 <- c(0.050, 0.005)
X2 <- c(0, 0.005)
DF7 <- cbind(X1, X2)
propagate(expr = EXPR7, data = DF7, type = "stat",
          nsim = 10000, check = FALSE)
## This will give exactly the same results as in 
## 9.4.2.2.7, Table 8, x1 = 0.050.

######### GUM 2008 Supplement 1 (4) #######################
## Example from 9.5 of the GUM 2008 Supplement 1
## (see 'References'), gauge block calibration.
## Assignment of PDF's as in Table 10 of 9.5.2.1.
EXPR8 <- expression(Ls + D + d1 + d2 - Ls * (da * (t0 + Delta) + as * dt) - Lnom)
Lnom <- 50000000
Ls <- propagate:::rst(10000, mean = 50000623, sd  = 25, df = 18)
D <- propagate:::rst(10000, mean = 215, sd = 6, df = 24)
d1 <- propagate:::rst(10000, mean = 0, sd = 4, df = 5)
d2 <- propagate:::rst(10000, mean = 0, sd = 7, df = 8)
as <- runif(10000, 9.5E-6, 13.5E-6)
t0 <- rnorm(10000, -0.1, 0.2)
Delta <- propagate:::rarcsin(10000, -0.5, 0.5)
da <- propagate:::rctrap(10000, -1E-6, 1E-6, 0.1E-6)
dt <- propagate:::rctrap(10000, -0.050, 0.050, 0.025)
DF8 <- cbind(Ls, D, d1, d2, as, t0, Delta, da, dt, Lnom)
propagate(expr = EXPR8, data = DF8, type = "sim", cov = "diag", alpha = 0.01,
          nsim = 10000, check = FALSE)
## This gives "fairly" the same results as in 9.5.4.2, Table 11,
## especially when comparing u.MC with u.2!
## GUM 2008 gives 32 and 36, respectively.
}   
}

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

