summary.propagate <- function(object, normality = FALSE, ...)
{
  ## print error propagation results
  cat(red("Results from uncertainty propagation:\n"))
  print(object$prop)
  
  ## print simulation results
  if (is.finite(object$sim[1])) {
    cat(red("Results from Monte Carlo simulation:\n"))
    print(object$sim)
  }
  
  ## print WS degrees-of-freedom
  cat(red("Welch-Satterthwaite degrees of freedom:\n"))
  print(object$ws.df)
  
  ## print coverage factor
  cat(red("Coverage factor (k):\n"))
  print(object$k)
  
  ## print coverage factor
  cat(red("Expanded uncertainty:\n"))
  print(object$u.exp)
  
  ## print symbolic derivatives of gradient
  if (is.matrix(object$gradient)) {
    cat(red("Symbolic gradient matrix:\n"))
    print(as.character(object$gradient))
  }
  
  ## print evaluated derivatives of gradient
  if (is.matrix(object$gradient)) {
    cat(red("Evaluated gradient matrix (Sensitivities):\n"))
    print(object$evalGrad)
  }
   
  ## print rel. contribution
  cat(red("Relative contribution:\n"))
  print(signif(diag(object$rel.contr), 3))
  
  ## print symbolic derivatives of hessian
  if (is.matrix(object$hessian)) {
    cat(red("Symbolic Hessian matrix:\n"))
    NCOL <- ncol(object$hessian)    
    STR <- as.character(object$hessian) 
    SPLIT <- split(STR, as.factor(rep(1:NCOL, NCOL)))
    for (i in 1:length(SPLIT)) print(SPLIT[[i]])
  }
    
  ## print evaluated derivatives of hessian
  if (is.matrix(object$hessian)) {
    cat(red("Evaluated Hessian matrix:\n"))
    NCOL <- ncol(object$hessian)
    EVAL <- object$evalHess 
    SPLIT <- split(EVAL, as.factor(rep(1:NCOL, NCOL)))
    for (i in 1:length(SPLIT)) print(SPLIT[[i]]) 
  }
  
  ## print input covariance matrix
  cat(red("Covariance matrix (Input):\n"))
  print(signif(object$covMat, 3))
  
  ## print MC-derived covariance matrix
  cat(red("Covariance matrix (MC):\n"))
  print(signif(object$covMat.MC), 3)
  
  ## print Frobenius norm of covariance matrices
  DIM <- ncol(object$covMat)
  NSIM <- object$nsim
  frobCut <- 3 * sqrt(DIM * (DIM - 1) / (2 * NSIM)) 
  cat(red("Frobenius Norm of input vs. MC-derived covariance matrices, should be < ", signif(frobCut, 3), ":\n"))
  print(signif(object$frobCOV), 3)
  
  ## print input correlation matrix
  cat(red("Correlation matrix (Input):\n"))
  print(signif(object$corMat, 3))
  
  ## print MC-derived correlation matrix
  cat(red("Correlation matrix (MC):\n"))
  print(signif(object$corMat.MC), 3)
  
  ## print Frobenius norm of correlation matrices
  cat(red("Frobenius Norm of input vs. MC-derived correlation matrices, should be < ", signif(frobCut, 3), ":\n"))
  print(signif(object$frobCOR), 3)
  
  ## comparison input and MC-output t-distributions
  if (!is.null(object$check)) {
    cat(red("Comparison of Copula-input and MC-derived moments/DOFs from fitted t-distribution:\n"))
    OUT <- rbind(object$data[1, ], object$check[1, ], object$data[2, ], object$check[2, ], object$data[3, ], object$check[3, ])
    rownames(OUT) <- c("mean.input", "mean.MC", "u.input", "u.MC", "df.input", "df.MC")
    print(OUT)
  }
  
  ## Skewness and excess kurtosis of evaluated MC simulations
  cat(red("Skewness / Excess Kurtosis of MC evaluations:\n"))
  cat(skewness(object$resSIM), "/", kurtosis(object$resSIM), "\n")
  
  if (normality) {
    cat(red("Repeated subsampling tests for normality:\n"))
    ## Shapiro-Wilk test for normality of MC distribution
    PVALS <- replicate(2000, shapiro.test(sample(object$resSIM, 5000))$p.value)
    medianP <- exp(median(log(PVALS), na.rm = TRUE))
    cat("Shapiro-Wilk median log(p-value), 2000r x 5000s: ", medianP)
    if (medianP >= 0.05) cat(" => normal\n") else cat(" => non-normal\n")
 
    ## Anderson-Darling test for normality of MC distribution
    PVALS <- replicate(1000, ad.test(sample(object$resSIM, 10000))$p.value)
    medianP <- exp(median(log(PVALS), na.rm = TRUE))
    cat("Anderson-Darling median log(p-value), 1000r x 10000s: ", medianP)
    if (medianP >= 0.05) cat(" => normal\n") else cat(" => non-normal\n")
  
    ## Lilliefors test for normality of MC distribution
    PVALS <- replicate(1000, lillie.test(sample(object$resSIM, 10000))$p.value)
    medianP <- exp(median(log(PVALS), na.rm = TRUE))
    cat("Lilliefors median log(p-value), 1000r x 10000s: ", medianP)
    if (medianP >= 0.05) cat(" => normal\n") else cat(" => non-normal\n")
  }    
}