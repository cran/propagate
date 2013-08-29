########### generating random deviates from distributions #############
## random samples from the skew-normal distribution, taken from package 'VGAM'
rsn <- function(n, location = 0, scale = 1, shape = 0) 
{
  rho <- shape/sqrt(1 + shape^2)
  u0 <- rnorm(n)
  v <- rnorm(n)
  u1 <- rho * u0 + sqrt(1 - rho^2) * v
  res <- location + scale * ifelse(u0 >= 0, u1, -u1)
  res[scale <= 0] <- NA
  res
}

## random samples from the generalized normal distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rgnorm <- function(n, location = 0, scale = 1, shape = 0) 
{
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
  res <- location + (scale * (1 - exp((sqrt(2) * shape * erf.inv(1 - 2 * runif(n))))))/shape
  res
}

## random samples from the scales and shifted t-distribution, 
## taken from the 'metRology' package
rst <- function(n, mean = 0, sd = 1, df = 2) 
{
  mean + sd * rt(n, df = df)
}

## random samples from the Gumbel distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rgumbel <- function(n, location = 0, scale = 1) 
{
  res <- location - scale * log(-log(runif(n)))
  res[scale <= 0] <- NaN
  res
}

## random samples from the Johnson SU distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rJSU <- function(n, xi = 0, lambda = 1, gamma = 1, delta = 1) 
{
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
  res <- xi + lambda * sinh((-gamma + sqrt(2) * erf.inv(-1 + 2 * runif(n)))/delta)
  res
}

## random samples from the Johnson SB distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rJSB <- function(n, xi = 0, lambda = 1, gamma = 1, delta = 1) 
{
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
  res <- xi + lambda/(1 + exp((gamma - sqrt(2) * erf.inv(-1 + 2 * runif(n)))/delta))
  res
}

## random samples from the three-parameter weibull distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rweibull2 <- function(n, location = 0, shape = 1, scale = 1) 
{
  res <- location + scale * (-log(1 - runif(n))^(1/shape))
  res
}

## random samples from the four-parameter beta distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rbeta2 <- function(n, alpha1 = 1, alpha2 = 1, a = 0, b = 0)  
{
  res <- a + qbeta(runif(n), shape1 = alpha1, shape2 = alpha2, lower.tail = TRUE, log.p = FALSE)
  res <- rescale(res, a, b)
  res
}

## random samples from the triangular distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rtriang <- function(n, a = 0, b = 0.5, c = 1)  
{
  res <- numeric(length(n))
  x <- runif(n)
  CRIT <- (b - a)/(c - a)
  res[x >= 0 & x <= CRIT] <- a + sqrt(((b - a) * (c - a)) * x[x >= 0 & x <= CRIT])
  res[CRIT < x & x <= 1] <- c - sqrt(((c - b) * (c - a)) * (1 - x[CRIT < x & x <= 1]))
  res
}

## random samples from the trapezoidal distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rtrap <- function(n, a = 0, b = 1, c = 2, d = 3)  
{
  res <- numeric(length(n))
  x <- runif(n)
  CRIT1 <- (a - b)/(a + b - c - d)
  CRIT2 <- (a + b - 2 * c)/(a + b - c - d)
  res[x > 0 & x < CRIT1] <- a + sqrt(((b - a) * (-a - b + c + d)) * x[x > 0 & x < CRIT1])
  res[CRIT1 <= x & x < CRIT2] <- 0.5 * (a + b + (- a - b + c + d) * x[CRIT1 <= x & x < CRIT2])
  res[CRIT2 <= x & x < 1] <- d - sqrt(((d - c) * (-a - b + c + d)) * (1 - x[CRIT2 <= x & x < 1]))
  res
}

## random samples from the Laplacian distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rlaplace <- function(n, mean = 0, sigma = 1)   
{
  RUNIF <- runif(n)
  res <- mean - (sigma * log(1 - (-1 + 2 * RUNIF) * sign(-1 + 2 * RUNIF))) * sign(-1 + 2 * RUNIF) 
  res
}

## random samples from the Arcsine distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rarcsin <- function(n, a = 2, b = 1)   
{
  res <- a + (b - a) * (sin(pi * runif(n)/2)^2)
  res 
}

## random samples from the von Mises distribution,  
## using "Rejection Sampling" and runif as the enveloping distribution
rmises <- function(n, mu = 1, kappa = 3)  
{
  nSAMPLE <- 0 
  N <- n
  OUT <- NULL
  
  while(nSAMPLE < n) {
    # (1) sample from rectangular distribution 
    u <- runif(N)
    # (2) Sample from enveloping rectangular distribution from mu - pi/mu + pi
    x <- runif(N, mu - pi, mu + pi)
    # (3) Accept candidate value, if probability u is smaller or equal than the density g(x)
    # at x divided by the density of f(x) * A (scale) 
    if (nSAMPLE == 0) {
      maxDUNIF <- max(dunif(x, mu - pi, mu + pi), na.rm = TRUE)
      maxDMISES <- max(dmises(x, mu = mu, kappa = kappa), na.rm = TRUE)
      A <- 1.1 * maxDMISES/maxDUNIF      
    }
    VAL <- 1/A * dmises(x, mu = mu, kappa = kappa)/dunif(x, mu - pi, mu + pi)
    ## increase OUT   
    OUT <- c(OUT, x[u <= VAL])
    ## calculate difference to 'n'
    DIFF <- n - length(OUT)
    ## end criterion for 'while'
    nSAMPLE <- length(OUT)
    ## if nSAMPLE < n, restart with twice the remaining values
    N <- 2 * DIFF   
  }  
  
  return(OUT[1:n])
}

## random samples from the curvilinear trapezoidal distribution,  
## using "Rejection Sampling" and runif as the enveloping distribution
rctrap <- function(n, a = 0, b = 1, d = 0.1)
{
  nSAMPLE <- 0 
  N <- n
  OUT <- NULL
  
  while(nSAMPLE < n) {
    # (1) sample from rectangular distribution 
    u <- runif(N)
    # (2) Sample from enveloping rectangular distribution from min/max
    x <- runif(N, a, b)
    # (3) Accept candidate value, if probability u is smaller or equal than the density g(x)
    # at x divided by the density of f(x) * A (scale) 
    if (nSAMPLE == 0) {
      maxDUNIF <- max(dunif(x, a, b), na.rm = TRUE)
      maxDCTRAP <- max(dctrap(x, a = a, b = b, d = d), na.rm = TRUE)
      A <- 1.1 * maxDCTRAP/maxDUNIF      
    }
    VAL <- 1/A * dctrap(x, a = a, b = b, d = d)/dunif(x, a, b)
    ## increase OUT   
    OUT <- c(OUT, x[u <= VAL])
    ## calculate difference to 'n'
    DIFF <- n - length(OUT)
    ## end criterion for 'while'
    nSAMPLE <- length(OUT)
    ## if nSAMPLE < n, restart with twice the remaining values
    N <- 2 * DIFF   
  }  
  
  return(OUT[1:n])
}

## random samples from the generalized trapezoidal distribution,  
## using "Rejection Sampling" and runif as the enveloping distribution
rgtrap <- function(n, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
                   n3 = 2, alpha = 1) 
{
  nSAMPLE <- 0 
  N <- n
  OUT <- NULL
  
  while(nSAMPLE < n) {
    # (1) sample from rectangular distribution 
    u <- runif(N)
    # (2) Sample from enveloping rectangular distribution from min/max
    x <- runif(N, min, max)
    # (3) Accept candidate value, if probability u is smaller or equal than the density g(x)
    # at x divided by the density of f(x) * A (scale) 
    if (nSAMPLE == 0) {
      maxDUNIF <- max(dunif(x, min, max), na.rm = TRUE)
      maxDGTRAP <- max(dgtrap(x, min = min, mode1 = mode1, mode2 = mode2, max = max,
                              n1 = n1, n3 = n3, alpha = alpha), na.rm = TRUE)
      A <- 1.1 * maxDGTRAP/maxDUNIF      
    }
    VAL <- 1/A * dgtrap(x, min = min, mode1 = mode1, mode2 = mode2, max = max,
                        n1 = n1, n3 = n3, alpha = alpha)/dunif(x, min, max)
    ## increase OUT   
    OUT <- c(OUT, x[u <= VAL])
    ## calculate difference to 'n'
    DIFF <- n - length(OUT)
    ## end criterion for 'while'
    nSAMPLE <- length(OUT)
    ## if nSAMPLE < n, restart with twice the remaining values
    N <- 2 * DIFF   
  }  
  
  return(OUT[1:n])
}