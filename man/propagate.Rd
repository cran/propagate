\name{propagate}
\alias{propagate}
\encoding{latin1}

\title{Propagation of uncertainty using higher-order Taylor expansion and Monte Carlo simulation}

\description{
A general function for the calculation of uncertainty propagation by first-/second-order Taylor expansion and Monte Carlo simulation including covariances. Input data can be any symbolic/numeric differentiable expression and data based on replicates, summaries (mean & s.d.) or sampled from a distribution. Uncertainty propagation is based completely on matrix calculus accounting for full covariance structure. Monte Carlo simulation is conducted using multivariate normal or t-distributions with covariance structure.
}

\usage{
propagate(expr, data, type = c("stat", "raw", "sim"), second.order = TRUE,
          method = c("sym", "num"), do.sim = TRUE, dist.sim = c("norm", "t"), 
          df.t = NULL, use.cov = TRUE, nsim = 100000, use.eval = TRUE,
          alpha = 0.05, ...)  
}

\arguments{
  \item{expr}{an expression, such as \code{expression(x/y)}.}
  \item{data}{a dataframe or matrix containing either a) the replicates in columns when using \code{type = "raw"}, b) the means in first and standard deviations in second rows when using \code{type = "stat"} or c) sampled data generated from any of \R's \code{\link{distributions}} or those implemented in this package (\code{\link{rDistr}}), if \code{type = "sim"}. Column names must match the variable names.}
  \item{type}{either \code{"stat"} if means and standard deviations are supplied, \code{"raw"} if raw replicates are given or \code{"sim"} in case of simulated data.}
  \item{second.order}{logical. If \code{TRUE}, error propagation will be calculated with first- and second-order Taylor expansion. See 'Details'.}
  \item{method}{first- and second-order partial derivatives are either calculated \code{"sym"}bollically or \code{"num"}erically.}
  \item{do.sim}{logical. Should Monte Carlo simulation be applied?}
  \item{dist.sim}{\code{"norm"} will use a multivariate normal distribution for Monte Carlo simulation, \code{"t"} a multivariate t-distribution. See 'Details'.}
  \item{df.t}{degrees of freedom when using \code{dist.sim = "t"} and no raw observations are supplied (i.e. \code{type = "stat"} or \code{type = "sim"}). See 'Details'.}
  \item{use.cov}{logical or variance-covariance matrix with the same column descriptions as \code{data}. See 'Details'.}
  \item{nsim}{the number of Monte Carlo simulations to be performed, minimum is 10000.}  
  \item{use.eval}{logical. If \code{TRUE}, the Monte Carlo simulation is evaluated in an environment ("vectorized"), otherwise evaluation takes place "row-wise" (using \code{\link{apply}}). See 'Details'.}
  \item{alpha}{the 1 - confidence level.}
  \item{...}{other parameters to be supplied to future methods.}
}

\details{
The implemented methods are:\cr\cr
1) \bold{Monte Carlo simulation:}\cr
For each variable \eqn{\boldsymbol{m}} in \code{data}, simulated data \eqn{\boldsymbol{x = [X_1, X_2, \ldots, X_n]}} with \eqn{\boldsymbol{n}} = \code{nsim} samples is generated from a multivariate normal distribution \eqn{\boldsymbol{x_{m, n} \sim \mathcal{N}(\mu, \Sigma)}} or multivariate t-distribution \eqn{\boldsymbol{x_{m, n} \sim t(\mu, \Sigma, \nu)}} using means \eqn{\boldsymbol{\mu_m}} and covariance matrix \eqn{\boldsymbol{\Sigma}} constructed from the standard deviations \eqn{\boldsymbol{\sigma_m}} of each variable. All data is coerced into a new dataframe that has the same covariance structure as the initial \code{data}: \eqn{\boldsymbol{\Sigma(\mathtt{data}) = \Sigma(x_{m, n})}}. Each row \eqn{\boldsymbol{i = 1, \ldots, n}} of the simulated dataset \eqn{\boldsymbol{x_{m, n}}} is evaluated with \code{expr}, \eqn{\boldsymbol{y_i = f(x_{m, i})}}, and summary statistics (mean, sd, median, mad, confidence interval based on \code{alpha}) are calculated on \eqn{\boldsymbol{y}}. If sample sizes for raw replicated data when setting \code{type = "raw"} are small, simulation from a multivariate t-distribution is advocated (Possolo, 2010). If this is done, the resulting distribution can have very heavy tails resulting in less stringent estimates of  \eqn{\boldsymbol{\sigma_y}} and larger confidence intervals. In this case, the degrees of freedom used are \eqn{n-1}. However, when setting \code{type = "stat"} or \code{type = "sim"}, the degrees of freedom have to be supplied by the user in \code{df.t}.

2) \bold{Error propagation:}\cr
The propagated error is calculated by first-/second-order Taylor expansion accounting for full covariance structure using matrix algebra.\cr
The following transformations based on two variables \eqn{x_1, x_2} illustrate the equivalence of the matrix-based approach with well-known classical notations:\cr
\bold{First-order mean:} \eqn{\rm{E[y]} = f(\bar{x}_i)}\cr
\bold{First-order variance:} \eqn{\sigma_y^2 = { \color{red} \nabla_x\mathbf{C}_x\nabla_x^T}}:\cr
\deqn{{ \color{red}[\rm{j_1}\; \rm{j_2}] \left[ \begin{array}{cc} \sigma_1^2 & \sigma_1\sigma_2 \\ \sigma_2\sigma_1 & \sigma_2^2 \end{array} \right] \left[ \begin{array}{c} \rm{j_1} \\ \rm{j_2} \end{array} \right]} = \rm{j_1}^2 \sigma_1^2 + \rm{2 j_1 j_2} \sigma_1 \sigma_2 + \rm{j_2}^2 \sigma_2^2}
\deqn{= \underbrace{\sum_{i=1}^2 \rm{j_i}^2 \sigma_i^2 + 2\sum_{i=1\atop i \neq k}^2\sum_{k=1\atop k \neq i}^2 \rm{j_i j_k} \sigma_{ik}}_{\rm{classical\;notation}} = \frac{1}{1!} \left(\sum_{i=1}^2 \frac{\partial f}{\partial x_i} \sigma_i \right)^2}
\cr
\bold{Second-order mean:} \eqn{\rm{E}[y] = f(\bar{x}_i) + {\color{blue} \frac{1}{2}\rm{tr}(\mathbf{H}_{xx}\mathbf{C}_x)}}:\cr
\deqn{{ \color{blue} \frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \sigma_1^2 & \sigma_1\sigma_2 \\ \sigma_2\sigma_1 & \sigma_2^2 \end{array} \right]} = \frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} \sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 & \rm{h_1}\sigma_1\sigma_2 + \rm{h_2}\sigma_2^2 \\ \rm{h_3} \sigma_1^2 + \rm{h_4} \sigma_1\sigma_2 & \rm{h_3} \sigma_1\sigma_2 + \rm{h_4} \sigma_2^2 \end{array} \right]}
\deqn{ = \frac{1}{2}(\rm{h_1}\sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 + \rm{h_3}\sigma_1\sigma_2 + \rm{h_4}\sigma_2^2) = \frac{1}{2!} \left(\sum_{i=1}^2 \frac{\partial}{\partial x_i} \sigma_i \right)^2 \it f} 
\cr
\bold{Second-order variance:} \eqn{\sigma_y^2 = {\color{red} \nabla_x\mathbf{C}_x\nabla_x^T} + {\color{blue} \frac{1}{2}\rm{tr}(\mathbf{H}_{xx}\mathbf{C}_x\mathbf{H}_{xx}\mathbf{C}_x)]}}:\cr
\deqn{{\color{blue}\frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \rm{\sigma_1^2} & \rm{\sigma_1\sigma_2} \\ \rm{\sigma_2\sigma_1} & \rm{\sigma_2^2} \end{array} \right] \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \rm{\sigma_1^2} & \rm{\sigma_1\sigma_2} \\ \rm{\sigma_2\sigma_1} & \rm{\sigma_2^2} \end{array} \right]} = \ldots}
\deqn{= \frac{1}{2} (\rm{h_1}^2\sigma_1^4 + \rm{2h_1h_2}\sigma_1^3\sigma_2 + \rm{2h_1h_3}\sigma_1^3\sigma_2 + \rm{h_2}^2\sigma_1^2\sigma_2^2 + \rm{2h_2h_3}\sigma_1^2\sigma_2^2 + \rm{h_3}^2\sigma_1^2\sigma_2^2 + \rm{2h_1h_4}\sigma_1^2\sigma_2^2}
\deqn{+ \rm{2h_2h_4}\sigma_1\sigma_2^3 + \rm{2h_3h_4}\sigma_1\sigma_2^3 + \rm{h_4}^2\sigma_2^4 = \frac{1}{2} (\rm{h_1}\sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 + \rm{h_3}\sigma_1\sigma_2 + \rm{h_4}\sigma_2^2)^2}
\deqn{= \frac{1}{2!} \left( \left(\sum_{i=1}^2 \frac{\partial}{\partial x_i} \sigma_i \right)^2 \it f \right)^2}
\cr
with \eqn{\mathrm{E}(y)} = expectation of \eqn{y}, \eqn{\mathbf{\sigma_y^2}} = variance of \eqn{y}, \eqn{{\color{red} \nabla_x}} = the p x n gradient matrix with all partial first derivatives \eqn{{\color{red} \rm{j_i}}}, \eqn{\mathbf{C}_x} = the p x p covariance matrix, \eqn{{\color{blue}\mathbf{H}_{xx}}} the Hessian matrix with all partial second derivatives \eqn{{\color{blue} \rm{h_i}}}, \eqn{\sigma_i} = the uncertainties and \eqn{\rm{tr}(\cdot)} = the trace (sum of diagonal) of a matrix. Note that because the hessian matrices are symmetric matrices, \eqn{{\color{blue} \rm{h_2}} = {\color{blue} \rm{h_3}}}. For a detailed derivation, see 'References'.\cr
The second-order Taylor expansion corrects for bias in nonlinear expressions as the first-order Taylor expansion assumes linearity around \eqn{\bar{x}_i}.\cr 
Depending on the input expression, the error propagation may result in an error that is not normally distributed. The Monte Carlo simulation, starting with normal distributions of the variables, can clarify this. For instance, a high tendency from deviation of normality is encountered in formulas in which the error of the denominator is relatively large or in exponential models with a large error in the exponent.\cr 
For setups in which there is no symbolic derivation possible (i.e. \code{e <- expression(abs(x))} => "Function 'abs' is not in the derivatives table") one can switch to numeric differentiation by \code{method = "num"}.\cr
Setting \code{use.eval = TRUE} has significant speed enhancement in evaluating the Monte Carlo simulations (~ 10-fold),  much to the author's surprise!
}

\value{
A list with the following components:   
\item{datSIM}{a vector containing the \code{nsim} simulated multivariate values for each variable in column format.}  
\item{resSIM}{a vector containing the \code{nsim} values obtained from the row-wise expression evaluations \eqn{\boldsymbol{ f(x_{m, i})}} of the simulated data in \code{datSIM}.}  
\item{datPROP}{\code{nsim} values generated from a normal distribution with \eqn{\mu} and \eqn{\sigma} as calculated from the propagated error.}   
\item{gradient}{the symbolic gradient vector \eqn{\nabla_x} of partial first-order derivatives.} 
\item{evalGrad}{the evaluated gradient vector \eqn{\nabla_x} of partial first-order derivatives.}
\item{covMat}{the covariance matrix \eqn{\mathbf{C}_x} used for Monte Carlo simulation and error propagation.} 
\item{hessian}{the symbolic hessian matrix \eqn{\mathbf{H}_{xx}} of partial second-order derivatives.}
\item{evalHess}{the evaluated hessian matrix \eqn{\mathbf{H}_{xx}} of partial second-order derivatives.}
\item{prop}{a summary vector containing first-/second-order expectations and uncertainties as well as the confidence interval based on \code{alpha} from \code{datPROP}.} 
\item{sim}{a summary vector containing the mean, standard deviation, median, MAD as well as the confidence interval based on \code{alpha} from \code{datSIM}.} 
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\bold{Error propagation (in general):}\cr
An Introduction to error analysis.\cr
Taylor JR.\cr
University Science Books (1996), New York.

Evaluation of measurement data - Guide to the expression of uncertainty in measurement.\cr
JCGM 100:2008 (GUM 1995 with minor corrections).\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf}.

Evaluation of measurement data - Supplement 1 to the Guide to the expression of uncertainty in measurement - Propagation of distributions using a Monte Carlo Method.\cr
JCGM 101:2008.\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_101_2008_E.pdf}.

\bold{Higher-order Taylor expansion:}\cr
On higher-order corrections for propagating uncertainties.\cr
Wang CM & Iyer HK.\cr
\emph{Metrologia} (2005), \bold{42}: 406-410.

Propagation of uncertainty: Expressions of second and third order uncertainty with third and fourth moments.\cr
Mekid S & Vaja D.\cr
\emph{Measurement} (2008), \bold{41}: 600-609.

\bold{Matrix algebra for error propagation:}\cr
An Introduction to Error Propagation: Derivation, Meaning and Examples of Equation Cy = FxCxFx^t.\cr
\url{www.nada.kth.se/~kai-a/papers/arrasTR-9801-R3.pdf}.

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
## From summary data: 
EXPR1 <- expression(x/y)
x <- c(5, 0.01)
y <- c(1, 0.01)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
RES1$prop
RES1$sim

## From raw data:
EXPR2 <- expression(x/y)
x <- c(2, 2.1, 2.2, 2, 2.3, 2.1)
y <- c(4, 4, 3.8, 4.1, 3.1, NA)
DF2 <- cbind(x, y)  
RES2 <- propagate(expr = EXPR2, data = DF2, type = "raw",  
                  do.sim = TRUE, verbose = TRUE)
RES2$prop
RES2$sim

## Compare to using a multivariate t-distribution
## because of low sample size => larger confidence
## intervals.
RES2b <- propagate(expr = EXPR2, data = DF2, type = "raw",  
                  do.sim = TRUE, verbose = TRUE, dist.sim = "t")
RES2$sim
RES2b$sim

################# GUM 2008 (1) ########################
## Example in Annex H.1 from the GUM 2008 manual
## (see 'References'), an end gauge calibration
## study. At first, we will only use the first-order 
## Taylor expansion.
EXPR3 <- expression(ls + d - ls * (da * the + as * dt))
ls <- c(50000623, 25)
d <- c(215, 9.7)
da <- c(0, 0.58E-6)
the <- c(-0.1, 0.41)
as <- c(11.5E-6, 1.2E-6)
dt <- c(0, 0.029)
DF3 <- cbind(ls, d, da, the, as, dt)
RES3 <- propagate(expr = EXPR3, data = DF3, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, 
                  second.order = FALSE)
RES3$prop
RES3$sim
## propagate: sd.1 = 31.71 
## GUM H.1.4/H.6c: u = 32  

## Proof that covariance of Monte-Carlo
## simulated dataset is the same as from 
## initial data.
RES3$covMat
cov(RES3$datSIM)
all.equal(RES3$covMat, cov(RES3$datSIM))

## Expanded uncertainty GUM H.1.6
## with 16 degrees of freedom.
qt(0.005, 16, lower.tail = FALSE) * 31.71
## propagate: 92.62
## GUM H.1.6: 93  

## Second-order terms GUM H.1.7.
RES4 <- propagate(expr = EXPR3, data = DF3, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, 
                  second.order = TRUE)
RES4$prop
RES4$sim
## propagate: sd.2 = 33.91115
## GUM H.1.7: u = 34.
## Also similar to the non-matrix-based approach
## in Wang et al. (2005, page 408): u1 = 33.91115.
## NOTE: After second-order correction, uncertainty is more
## similar to the value obtained from Monte Carlo simulation!

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

## H.2.3 Approach 1 using mean values and
## standard uncertainties:
EXPR6a <- expression((V/I) *  cos(phi)) ## R
EXPR6b <- expression((V/I) *  sin(phi)) ## X
EXPR6c <- expression(V/I) ## Z
MEAN6 <- colMeans(H.2)
SD6 <- sqrt(colVarsC(H.2))
DF6 <- rbind(MEAN6, SD6)
COV6ab <- cov(H.2) ## covariance matrix of V, I, phi
COV6c <- cov(H.2[, 1:2])  ## covariance matrix of V, I

RES6a <- propagate(expr = EXPR6a, data = DF6, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, use.cov = COV6ab,
                  second.order = TRUE)

RES6b <- propagate(expr = EXPR6b, data = DF6, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, use.cov = COV6ab,
                  second.order = TRUE)
                  
RES6c <- propagate(expr = EXPR6c, data = DF6[, 1:2], type = "stat", 
                  do.sim = TRUE, verbose = TRUE, use.cov = COV6c,
                  second.order = TRUE)

## This gives exactly the same values of mean and sd/sqrt(5)
## as given in Table H.4.
RES6a$prop
RES6b$prop
RES6c$prop

#################### GUM 2008 (3) #######################
## Example in Annex H.3 from the GUM 2008 manual
## (see 'References'), calibration of a thermometer.
## For this example, we can use the predict.lm function
## of R directly to return the error of predicted values.
data(H.3)
LM <- lm(bk ~ I(tk - 20), data = H.3)
## This will give exactly the same values as in H.3.3.
summary(LM)

## This will give exactly the same values as the
## fourth column ("Predicted correction") of Table H.6.
predict(LM)

## This will give exactly the same values as the
## fifth column ("Difference...") of Table H.6.
H.3$bk - predict(LM)

## Uncertainty in a predicted value. This will give
## exactly the values in H.3.4.
predict(LM, newdata = data.frame(tk = 30), se.fit = TRUE)

######### GUM 2008 Supplement 1 (1) #######################
## Example from 9.2.2 of the GUM 2008 Supplement 1
## (see 'References'), normally distributed input
## quantities. Assign values as in 9.2.2.1.
EXPR7 <- expression(X1 + X2 + X3 + X4)
X1 <- c(0, 1)
X2 <- c(0, 1)
X3 <- c(0, 1)
X4 <- c(0, 1)
DF7 <- cbind(X1, X2, X3, X4)
RES7 <- propagate(expr = EXPR7, data = DF7, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, nsim = 1E5)
## This will give exactly the same results as in 
## 9.2.2.6, Table 2.
RES7$prop
RES7$sim

######### GUM 2008 Supplement 1 (2) #######################
## Example from 9.3 of the GUM 2008 Supplement 1
## (see 'References'), mass calibration.
## Formula 24 in 9.3.1.3 and values as in 9.3.1.4, Table 5.
EXPR8 <- expression((Mrc + dMrc) * (1 + (Pa - Pa0) * ((1/Pw) - (1/Pr))) - Mnom)
Mrc <- rnorm(1E5, 100000, 0.050)
dMrc <- rnorm(1E5, 1.234, 0.020)
Pa <- runif(1E5, 1.10, 1.30)  ## E(Pa) = 1.2, (b-a)/2 = 0.1 
Pw <- runif(1E5, 7000, 9000)  ## E(Pw) = 8000, (b-a)/2 = 1000
Pr <- runif(1E5, 7950, 8050) ## E(Pr) = 8000, (b-a)/2 = 50
Pa0 <- 1.2 
Mnom <- 100000
DF8 <- cbind(Mrc, dMrc, Pa, Pw, Pr, Pa0, Mnom)
RES8 <- propagate(expr = EXPR8, data = DF8, type = "sim",
                  do.sim = TRUE, verbose = TRUE, nsim = 1E5)
## This will give exactly the same results as in 
## 9.3.2.3, Table 6
RES8$prop
RES8$sim
 
######### GUM 2008 Supplement 1 (3) #######################
## Example from 9.4 of the GUM 2008 Supplement 1
## (see 'References'), comparioson loss in microwave
## power meter calibration, zero covariance.
## Formula 28 in 9.4.1.5 and values as in 9.4.1.7.
EXPR9 <- expression(X1^2 - X2^2)
X1 <- c(0.050, 0.005)
X2 <- c(0, 0.005)
DF9 <- cbind(X1, X2)
RES9a <- propagate(expr = EXPR9, data = DF9, type = "stat",
                  do.sim = TRUE, verbose = TRUE, nsim = 1E5)
## This will give exactly the same results as in 
## 9.4.2.2.7, Table 8, x1 = 0.050.
RES9a$prop
RES9a$sim

## Using covariance matrix with r(x1, x2) = 0.9
## We convert to covariances using cor2cov.
COR9 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
COV9 <- cor2cov(COR9, c(0.005^2, 0.005^2))
colnames(COV9) <- c("X1", "X2")
rownames(COV9) <- c("X1", "X2")
RES9b <- propagate(expr = EXPR9, data = DF9, type = "stat",
                   use.cov = COV9, do.sim = TRUE, 
                   verbose = TRUE, nsim = 1E5)
## This will give exactly the same results as in 
## 9.4.3.2.1, Table 9, x1 = 0.050.
RES9b$prop
RES9b$sim

######### GUM 2008 Supplement 1 (4) #######################
## Example from 9.5 of the GUM 2008 Supplement 1
## (see 'References'), gauge block calibration.
## Assignment of PDF's as in Table 10 of 9.5.2.1.
EXPR10 <- expression(Ls + D + d1 + d2 - Ls *(da *(t0 + Delta) + as * dt) - Lnom)
Lnom <- 50000000
Ls <- propagate:::rst(100000, mean = 50000623, sd  = 25, df = 18)
D <- propagate:::rst(100000, mean = 215, sd = 6, df = 25)
d1 <- propagate:::rst(100000, mean = 0, sd = 4, df = 5)
d2 <- propagate:::rst(100000, mean = 0, sd = 7, df = 8)
as <- runif(100000, 9.5E-6, 13.5E-6)
t0 <- rnorm(100000, -0.1, 0.2)
Delta <- propagate:::rarcsin(100000, -0.5, 0.5)
da <- propagate:::rctrap(100000, -1E-6, 1E-6, 0.1E-6)
dt <- propagate:::rctrap(100000, -0.050, 0.050, 0.025)
DF10 <- cbind(Ls, D, d1, d2, as, t0, Delta, da, dt, Lnom)
RES10 <- propagate(expr = EXPR10, data = DF10, type = "sim",
                   use.cov = FALSE, verbose = TRUE, 
                   alpha = 0.01)
## This gives the same results as in 9.5.4.2, Table 11.
## However: results are exacter than in the GUM 2008
## manual, especially when comparing sd(Monte Carlo)
## = 34.49 and sd(second-order Taylor expansion) 
## = 34.45, while GUM 2008 gives 32 and 36, respectively.
RES10$prop
RES10$sim
}   

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

