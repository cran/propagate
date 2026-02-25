sobol <- function(prop = NULL, A = NULL, B = NULL, expr = NULL, method = c("jansen", "sobol", "saltelli", "homma")) {
  
  method <- match.arg(method)
  
  ## external A & B & expression for gaugun, e.g. Ishigami
  if (is.matrix(A)) {
    if (!is.matrix(B) || !is.expression(expr)) stop("'A','B' must be matrices & 'expr' be an expression")
    EXPR <- expr
    DATA <- A
  
  } else {
    ## extract from propagate object
    if (!inherits(prop, "propagate")) stop("Input must be a 'propagate' object!")
    
    EXPR <- prop$expr
    DATA <- prop$data
    
    A <- prop$datSIM
    B <- prop$datSIM2
  }
  
  d <- ncol(A)
  
  ## hybrid matrices
  A_B <- vector("list", d)
  for (i in 1:d) {
    tmp <- A
    tmp[, i] <- B[, i]
    A_B[[i]] <- tmp
  }
  
  ## model evaluation
  Y_A  <- propagate(EXPR, A, type = "sim", cov = "full", check = FALSE, exp.uc = "1st")$resSIM
  Y_B  <- propagate(EXPR, B, type = "sim", cov = "full", check = FALSE, exp.uc = "1st")$resSIM
  Y_AB <- lapply(A_B, function(x) propagate(EXPR, x, type = "sim", cov = "full", check = FALSE, exp.uc = "1st")$resSIM)
  
  # Sobol estimation
  VY <- var(c(Y_A, Y_B))
  
  S  <- numeric(d)
  ST <- numeric(d)
  
  for (i in 1:d) {
    if (method == "sobol") {
      S[i]  <- mean(Y_B * (Y_AB[[i]] - Y_A)) / VY
      ST[i] <- NA
    }
    
    if (method == "saltelli") {
      S[i]  <- mean(Y_B * (Y_AB[[i]] - Y_A)) / VY
      ST[i] <- mean((Y_A - Y_AB[[i]])^2) / (2 * VY)
    }
    
    if (method == "homma") {
      S[i]  <- NA
      ST[i] <- mean(Y_A * (Y_A - Y_AB[[i]])) / VY
    }
    
    if (method == "jansen") {
      S[i]  <- 1 - mean((Y_B - Y_AB[[i]])^2) / (2 * VY)
      ST[i] <- mean((Y_A - Y_AB[[i]])^2) / (2 * VY)
    }
  }
  
  names(S)  <- names(ST) <- colnames(DATA)
  return(list(S = S, ST = ST))
}
