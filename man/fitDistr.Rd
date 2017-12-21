\name{fitDistr}
\alias{fitDistr}

\title{Fitting distributions to observations/Monte Carlo simulations}

\description{
This function fits 29 different continuous distributions by (weighted) NLS to the histogram of Monte Carlo simulation results as obtained by \code{\link{propagate}} or any other vector containing large-scale observations. Finally, the fits are sorted by ascending \code{\link{AIC}}. 
}

\usage{
fitDistr(object, nbin = 100, weights = NULL, verbose = TRUE, plot = TRUE,  ...)
}

\arguments{
  \item{object}{an object of class 'propagate' or a vector containing observations.}
  \item{nbin}{in case of \code{type = "hist"}, the number of bins in the histogram.}
  \item{weights}{numeric or logical. Either a numeric vector of weights, or if \code{TRUE}, the distributions are fitted with weights = 1/(counts per bin).}
  \item{verbose}{logical. If \code{TRUE}, steps of the analysis are printed to the console.}
  \item{plot}{logical. If \code{TRUE}, a plot with the "best" distribution (in terms of lowest AIC) on top of the histogram or kernel density curve is displayed.}
  \item{...}{other parameters to be passed to \code{\link{density}} or \code{\link{histogram}} .}  
}

\details{
Fits the following 29 distributions using (weighted) residual sum-of-squares as the minimization criterion for \code{\link{nls.lm}}:\cr
1) Normal distribution (\code{dnorm}) => \url{http://en.wikipedia.org/wiki/Normal_distribution}\cr
2) Skewed-normal distribution (\code{propagate:::dsn}) => \url{http://en.wikipedia.org/wiki/Skew_normal_distribution}\cr
3) Generalized normal distribution (\code{propagate:::dgnorm}) => \url{http://en.wikipedia.org/wiki/Generalized_normal_distribution}\cr 
4) Log-normal distribution (\code{dlnorm}) => \url{http://en.wikipedia.org/wiki/Log-normal_distribution}\cr 
5) Scaled and shifted t-distribution (\code{propagate:::dst}) => GUM 2008, Chapter 6.4.9.2.\cr
6) Logistic distribution (\code{dlogis}) => \url{http://en.wikipedia.org/wiki/Logistic_distribution}\cr 
7) Uniform distribution (\code{dunif}) => \url{http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)}\cr 
8) Triangular distribution (\code{propagate:::dtriang}) => \url{http://en.wikipedia.org/wiki/Triangular_distribution}\cr 
9) Trapezoidal distribution (\code{propagate:::dtrap}) => \url{http://link.springer.com/article/10.1007\%2Fs001840200230}\cr 
10) Curvilinear Trapezoidal distribution (\code{propagate:::dctrap}) => GUM 2008, Chapter 6.4.3.1\cr
11) Generalized Trapezoidal distribution (\code{propagate:::dgtrap}) => \url{http://www.seas.gwu.edu/~dorpjr/Publications/JournalPapers/Metrika2003VanDorp.pdf}\cr
12) Gamma distribution (\code{dgamma}) => \url{http://en.wikipedia.org/wiki/Gamma_distribution}\cr 
13) Inverse Gamma distribution (\code{propagate:::dinvgamma}) => \url{http://en.wikipedia.org/wiki/Inverse-gamma_distribution}\cr 
14) Cauchy distribution (\code{dcauchy}) => \url{http://en.wikipedia.org/wiki/Cauchy_distribution}\cr 
15) Laplace distribution (\code{propagate:::dlaplace}) => \url{http://en.wikipedia.org/wiki/Laplace_distribution}\cr 
16) Gumbel distribution (\code{propagate:::dgumbel}) => \url{http://en.wikipedia.org/wiki/Gumbel_distribution}\cr 
17) Johnson SU distribution (\code{propagate:::dJSU}) => \url{http://en.wikipedia.org/wiki/Johnson_SU_distribution}\cr 
18) Johnson SB distribution (\code{propagate:::dJSB}) => \url{http://www.mathwave.com/articles/johnson_sb_distribution.html}\cr 
19) Three-parameter Weibull distribution (\code{propagate:::dweibull2}) => \url{http://en.wikipedia.org/wiki/Weibull_distribution}\cr 
20) Two-parameter beta distribution (\code{dbeta2}) => \url{https://en.wikipedia.org/wiki/Beta_distribution#Two_parameters_2}\cr
21) Four-parameter beta distribution (\code{propagate:::dbeta2}) => \url{https://en.wikipedia.org/wiki/Beta_distribution#Four_parameters_2}\cr 
22) Arcsine distribution (\code{propagate:::darcsin}) => \url{http://en.wikipedia.org/wiki/Arcsine_distribution}\cr 
23) Von Mises distribution (\code{propagate:::dmises}) => \url{http://en.wikipedia.org/wiki/Von_Mises_distribution}\cr 
24) Inverse Gaussian distribution (\code{propagate:::dinvgauss}) => \url{http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}\cr 
25) Generalized Extreme Value distribution (\code{propagate:::dgevd}) => \url{http://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}\cr 
26) Rayleigh distribution (\code{propagate:::drayleigh}) => \url{http://en.wikipedia.org/wiki/Rayleigh_distribution}\cr
27) Chi-square distribution (\code{dchisq}) => \url{https://en.wikipedia.org/wiki/Chi-squared_distribution}\cr
28) Exponential distribution (\code{dexp}) => \url{https://en.wikipedia.org/wiki/Exponential_distribution}\cr
28) F-distribution (\code{df}) => \url{https://en.wikipedia.org/wiki/F-distribution}\cr

Distributions 3), 17), 18), 19) and 21) are sometimes hard to fit because the start parameters are not readily deducible from the kernel density estimates or some parameters are highly sensitive to shape changes. For these five cases, a grid of starting values with different magnitudes is used to obtain the best parameter combination with respect to lowest residual sum-of-squares ("brute force" approach).

The goodness-of-fit (GOF) is calculated with \code{\link{AIC}} from the (weighted) log-likelihood of the fit:
\deqn{\rm{ln}(L) = 0.5 \cdot \left(-N \cdot \left(\rm{ln}(2\pi) + 1 + \rm{ln}(N) - \sum_{i=1}^n log(w_i) + \rm{ln}\left(\sum_{i=1}^n w_i \cdot x_i^2\right) \right) \right) \; , \; AIC = 2k - 2\rm{ln}(L)}
with \eqn{x_i} = the residuals from the NLS fit, \eqn{N} = the length of the residual vector, \eqn{k} = the number of parameters of the fitted model and \eqn{w_i} = the weights.

In contrast to some other distribution fitting softwares (i.e. Easyfit, Mathwave) that use residual sum-of-squares/Anderson-Darling/Kolmogorov-Smirnov statistics as GOF measures, the application of AIC accounts for increasing number of parameters in the distribution fit and therefore compensates for overfitting. Hence, this approach is more similar to ModelRisk (Vose Software) and as employed in \code{fitdistr} of the 'MASS' package.
Another application is to identify a possible distribution for the raw data prior to using Monte Carlo simulations from this distribution. However, a decent number of observations should be at hand in order to obtain a realistic estimate of the proper distribution. See 'Examples'.\cr
The code for the density functions is in file "distr-densities.R".

IMPORTANT: It can be feasible to set \code{weights = TRUE} in order to give more weight to bins with low counts. See 'Examples'.

ALSO: Distribution fitting is highly sensitive to the number of defined histogram bins, so it is advisable to change this parameter and inspect if the order of fitted distributions remains stable!
}

\value{
A list with the following items:

\code{aic}: the ascendingly sorted AIC dataframe.\cr
\code{fit}: a list of the results from \code{\link{nls.lm}} for each distribution model.\cr
\code{bestfit}: the best model in terms of lowest AIC.\cr
\code{fitted}: the fitted values from the best model.\cr
\code{residuals}: the residuals from the best model.\cr
}
     
\author{
Andrej-Nikolai Spiess
}  

\references{
Continuous univariate distributions, Volume 1.\cr
Johnson NL, Kotz S and Balakrishnan N.\cr
\emph{Wiley Series in Probability and Statistics, 2.ed} (2004).

A guide on probability distributions.\cr
R-forge distributions core team.\cr
\url{http://dutangc.free.fr/pub/prob/probdistr-main.pdf}.

Univariate distribution relationships.\cr
Leemis LM and McQueston JT.\cr
\emph{The American Statistician} (2008), \bold{62}: 45-53.
}

\examples{
\dontrun{
## Linear example, small error
## => family of normal distributions.
EXPR1 <- expression(x + 2 * y)
x <- c(5, 0.01)
y <- c(1, 0.01)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
fitDistr(RES1)$aic

## Ratio example, larger error
## => family of slightly skewed distributions.
EXPR2 <- expression(x/2 * y)
x <- c(5, 0.1)
y <- c(1, 0.02)
DF2 <- cbind(x, y)
RES2 <- propagate(expr = EXPR2, data = DF2, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
fitDistr(RES2)$aic

## Exponential example, large error
## => family of stronger skewed distributions.
EXPR3 <- expression(x^(2 * y))
x <- c(5, 0.1)
y <- c(1, 0.1)
DF3 <- cbind(x, y)
RES3 <- propagate(expr = EXPR3, data = DF3, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
fitDistr(RES3)$aic

## Rectangular input distributions result
## in trapezoidal output distribution.
A <- runif(100000, 20, 25)
B <- runif(100000, 3, 3.5)
DF4 <- cbind(A, B)
EXPR4 <- expression(A + B)
RES4 <- propagate(EXPR4, data = DF4, type = "sim", 
                 use.cov = FALSE, do.sim = TRUE)
fitDistr(RES4)$aic        

## Fitting with 1/counts as weights.
EXPR5 <- expression(x + 2 * y)
x <- c(5, 0.05)
y <- c(1, 0.05)
DF5 <- cbind(x, y)
RES5 <- propagate(expr = EXPR5, data = DF5, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, weights = TRUE)
fitDistr(RES5)$aic
}
}   

\keyword{algebra}
\keyword{univariate}

