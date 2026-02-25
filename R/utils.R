########## visible ########################
makeGrad <- function(expr, order = NULL)
{
  VARS <- all.vars(expr)
  if (!is.null(order)) VARS <- VARS[order]
  FUN <- function(x) D(expr, x)
  vecGRAD <- sapply(VARS, FUN)
  vecGRAD <- matrix(vecGRAD, nrow = 1)  
  return(vecGRAD)  
} 

makeHess <- function(expr, order = NULL)
{
  VARS <- all.vars(expr)  
  if (!is.null(order)) VARS <- VARS[order]
  GRID <- expand.grid(VARS, VARS)    
  FUN <- function(x) D(D(expr, x[1]), x[2])
  vecHESS <- apply(GRID, 1, FUN)  
  matHESS <- matrix(vecHESS, ncol = length(VARS), byrow = TRUE)    
  return(matHESS)
} 

evalDerivs <- function(deriv, envir)
{
  if (missing(envir)) envir <- .GlobalEnv
  DIM <- dim(deriv)
  evalVEC <- sapply(deriv, eval, envir = envir)
  dim(evalVEC) <- DIM
  return(evalVEC)
}

kurtosis <- function (x, na.rm = FALSE) 
{
  if (is.matrix(x)) 
    apply(x, 2, kurtosis, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    n <- length(x)
    n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2) - 3
  }
  else if (is.data.frame(x)) 
    sapply(x, kurtosis, na.rm = na.rm)
  else kurtosis(as.vector(x), na.rm = na.rm)
}

skewness <- function (x, na.rm = FALSE) 
{
  if (is.matrix(x)) 
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x)) 
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

counter <- function (i) 
{
  if (i%%10 == 0) 
    cat(i)
  else cat(".")
  if (i%%50 == 0) 
    cat("\n")
  flush.console()
}

tr <- function(mat) sum(diag(mat), na.rm = TRUE)

rescale <- function (x, tomin, tomax) 
{
  if (missing(x) | missing(tomin) | missing(tomax)) {
    stop(paste("rescale: rescale(x, tomin, tomax)\n", "\twhere x is a numeric object and tomin and tomax\n is the range to rescale into", 
               sep = "", collapse = ""))
  }
  if (is.numeric(x) && is.numeric(tomin) && is.numeric(tomax)) {
    xrange <- range(x, na.rm = TRUE)
    if (xrange[1] == xrange[2]) 
      return(x)
    mfac <- (tomax - tomin)/(xrange[2] - xrange[1])
    return(tomin + (x - xrange[1]) * mfac)
  }
  else {
    warning("rescale: only numeric objects can be rescaled")
    return(x)
  }
}

print.interval <- function(x, ...)
{
  cat("[", x[1], ", ", x[2], "]", sep = "")
}

print.propagate <- function(x, ...)
{
  object <- x
  
  ## print error propagation results
  cat(red("Results from uncertainty propagation:\n"))
  print(object$prop)
  
  ## print simulation results
  cat(red("Results from Monte Carlo simulation:\n"))
  print(object$sim)
}

print.fitDistr <- function(x, ...)
{
  cat(red("Best fit is", names(x$fit)[[1]], "Distribution.\n"))
  cat("Parameters:\n")
  print(x$par[[1]])
  cat("Standard errors:\n")
  print(x$se[[1]])
  cat("Goodness of fit:\n")
  cat("BIC =", x$stat[1, "BIC"])
}

isFALSE <- function(x) is.logical(x) && length(x) == 1L && !is.na(x) && !x
isTRUE <- function(x) is.logical(x) && length(x) == 1L && !is.na(x) && x

ad.test <- function(x) 
{
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if (n < 8) 
    stop("sample size must be greater than 7")
  logp1 <- pnorm((x - mean(x))/sd(x), log.p = TRUE)
  logp2 <- pnorm(-(x - mean(x))/sd(x), log.p = TRUE)
  h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
  A <- -n - mean(h)
  AA <- (1 + 0.75/n + 2.25/n^2) * A
  if (AA < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
  }
  else if (AA < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
  }
  else if (AA < 0.6) {
    pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
  }
  else if (AA < 10) {
    pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
  }
  else pval <- 3.7e-24
  RVAL <- list(statistic = c(A = A), p.value = pval, method = "Anderson-Darling normality test", 
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

lillie.test <- function(x) 
{
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if (n < 5) 
    stop("sample size must be greater than 4")
  p <- pnorm((x - mean(x))/sd(x))
  Dplus <- max(seq(1:n)/n - p)
  Dminus <- max(p - (seq(1:n) - 1)/n)
  K <- max(Dplus, Dminus)
  if (n <= 100) {
    Kd <- K
    nd <- n
  }
  else {
    Kd <- K * ((n/100)^0.49)
    nd <- 100
  }
  pvalue <- exp(-7.01256 * Kd^2 * (nd + 2.78019) + 2.99587 * 
                  Kd * sqrt(nd + 2.78019) - 0.122119 + 0.974598/sqrt(nd) + 
                  1.67997/nd)
  if (pvalue > 0.1) {
    KK <- (sqrt(n) - 0.01 + 0.85/sqrt(n)) * K
    if (KK <= 0.302) {
      pvalue <- 1
    }
    else if (KK <= 0.5) {
      pvalue <- 2.76773 - 19.828315 * KK + 80.709644 * 
        KK^2 - 138.55152 * KK^3 + 81.218052 * KK^4
    }
    else if (KK <= 0.9) {
      pvalue <- -4.901232 + 40.662806 * KK - 97.490286 * 
        KK^2 + 94.029866 * KK^3 - 32.355711 * KK^4
    }
    else if (KK <= 1.31) {
      pvalue <- 6.198765 - 19.558097 * KK + 23.186922 * 
        KK^2 - 12.234627 * KK^3 + 2.423045 * KK^4
    }
    else {
      pvalue <- 0
    }
  }
  RVAL <- list(statistic = c(D = K), p.value = pvalue, method = "Lilliefors (Kolmogorov-Smirnov) normality test", 
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}