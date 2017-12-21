print.propagate <- function(x, ...)
{
  object <- x
  
  ## print error propagation results
  message("Results from uncertainty propagation:")
  print(object$prop)
  
  ## print simulation results
  if (length(x$resSIM) > 1) {
    message("Results from Monte Carlo simulation:")
    print(object$sim)
  }
}
