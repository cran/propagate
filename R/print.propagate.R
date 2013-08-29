print.propagate <- function(x, verbose = TRUE, ...)
{
  object <- x
  
  ## print error propagation results
  cat("Results from error propagation:\n")
  print(object$prop)
  
  ## print simulation results
  cat("\nResults from Monte Carlo simulation:\n")
  print(object$sim)
  
  if (verbose) {
    ## print covariance matrix
    cat("\nCovariance matrix:\n")
    print(object$covMat)
  
    ## print derivatives of gradient
    cat("\nSymbolic gradient matrix:\n")
    print(as.character(object$gradient))
  
    ## print derivatives of gradient
    cat("\nEvaluated gradient matrix:\n")
    print(object$evalGrad)
      
    ## print derivatives of hessian
    cat("\nSymbolic hessian matrix:\n")
    if (is.matrix(object$hessian)) {
      NCOL <- ncol(object$hessian)    
      STR <- as.character(object$hessian) 
      SPLIT <- split(STR, as.factor(rep(1:NCOL, NCOL)))
      for (i in 1:length(SPLIT)) print(SPLIT[[i]])
    } else cat("None available due to 'second.order = FALSE'.\n")
  
    ## print derivatives of hessian
    cat("\nEvaluated hessian matrix:\n")
    if (is.matrix(object$hessian)) {
      NCOL <- ncol(object$hessian)
      EVAL <- object$evalHess 
      SPLIT <- split(EVAL, as.factor(rep(1:NCOL, NCOL)))
      for (i in 1:length(SPLIT)) print(SPLIT[[i]])  
    } else cat("None available due to 'second.order = FALSE'.\n")
  }
}
