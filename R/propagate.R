propagate <- function(
expr, 
data, 
type = c("stat", "raw", "sim"),
cov = c("diag", "full"), 
df.tot = NULL,
k = NULL, 
nsim = 1000000,
alpha = 0.05,
exp.uc = c("1st", "2nd"),
check = TRUE,
...
)
{            
  op <- options(warn = -1)
  on.exit(options(op))
  type <- match.arg(type)
  if (is.matrix(cov)) cov <- cov else cov <- match.arg(cov)
  exp.uc <- match.arg(exp.uc)
  
  ## version 1.0-4: convert function to expression
  if (is.function(expr)) {
    ARGS <- as.list(args(expr))
    ARGS <- ARGS[-length(ARGS)]
    VARS <- names(ARGS)    
    expr <- body(expr)
    class(expr) <- "expression"
    isFun <- TRUE
  } else isFun <- FALSE
  
  ## check for correct expression and number of simulations
  if (!is.expression(expr)) stop("propagate: 'expr' must be an expression")
  if (nsim < 10000) stop("propagate: 'nsim' should be >= 10000 !")
  
  ## check for matching variable names
  if (!isFun) VARS <- all.vars(expr)  
  m <- match(VARS, colnames(data))
  
  ## version 1.0-8: if expr does not use all data columns, select
  if (length(m) != ncol(data)) data <- data[, m]
    
  if (any(is.na(m))) stop("propagate: variable names of input dataframe and expression do not match!")
  if (length(unique(m)) != length(m)) stop("propagate: some variable names are repetitive!")
  
  DATA <- data
  EXPR <- expr
  
  ## if covariance matrix is supplied, check for symmetry and matching names
  if (is.matrix(cov)) {
    if (NCOL(cov) != NROW(cov)) stop("propagate: 'cov' is not a symmetric matrix!")
    m <- match(colnames(cov), colnames(DATA))            
    if (any(is.na(m))) stop("propagate: names of input dataframe and var-cov matrix do not match!")             
    if (length(unique(m)) != length(m)) stop("propagate: some names of the var-cov matrix are repetitive!")        
  }
  
  ############## version 1.0-8: Create variables from input data, use the direct 'type's ##############
  
  ## 1) statistical summaries
  if (type == "stat") {
    meanVAL <- DATA[1, ]
    sdVAL <- DATA[2, ] 
    dfVAL <- if (nrow(DATA) == 3) DATA[3, ] else rep(1000, ncol(DATA))
    
    if (nrow(DATA) == 2 | nrow(DATA) == 3) cat("I see data with 2 or 3 rows ... type = 'stat' seems correct!\n")
    else cat("I don't see data with 2 or 3 rows ... type = 'stat' probably correct!\n")
    
    if (is.matrix(cov)) SIGMA <- cov 
    else if (cov == "diag") {
      SIGMA <- diag(sdVAL^2, nrow = length(VARS), ncol = length(VARS))    
      colnames(SIGMA) <- rownames(SIGMA) <- colnames(DATA)
    } 
    else if (cov == "full") stop("Cannot use full covariance matrix with summary data. Please supply an external one...")
    COR <- cov2cor(SIGMA)
  }
 
  ## 2) experimental samples
  if (type == "raw") {
    N <- apply(DATA, 2, function(x) sum(!is.na(x), na.rm = TRUE))
    meanVAL <- apply(DATA, 2, function(x) mean(x, na.rm = TRUE))
    sdVAL <- apply(DATA, 2, function(x) sd(x, na.rm = TRUE))/sqrt(N)
    dfVAL <- N - 1
    
    if (nrow(DATA) > 3 & nrow(DATA) < 1000) cat("I see data with 3 < rows < 1000 ... type = 'raw' seems correct!\n")
    else cat("I don't see data with 3 < rows < 1000 ... type = 'raw' probably correct!\n")
    
    if (is.matrix(cov)) {
      SIGMA <- cov
      COR <- cov2cor(SIGMA)
    }
    else if (cov == "diag") {
      SIGMA <- diag(sdVAL^2, nrow = length(VARS), ncol = length(VARS))    
      colnames(SIGMA) <- rownames(SIGMA) <- colnames(DATA)
      COR <- cov2cor(SIGMA)
    } 
    else if (cov == "full") {
      notNA <- !is.na(DATA)
      Npairs <- t(notNA) %*% notNA  
      SIGMA <- cov(DATA, use = "pairwise.complete.obs")/Npairs
      COR <- cor(DATA, use = "pairwise.complete.obs")
    }
  }
  
  ## 3) simulated distribution data
  if (type == "sim") {
    N <- nrow(DATA)
    meanVAL <- apply(DATA, 2, function(x) mean(x, na.rm = TRUE))
    sdVAL <- apply(DATA, 2, function(x) sd(x, na.rm = TRUE))
    dfVAL <- rep(N - 1, ncol(DATA))
    
    if (nrow(DATA) > 1000) cat("I see data with 1000 < rows ... type = 'sim' seems correct!\n")
    else cat("I don't see data with 1000 < rows ... type = 'sim' probably correct!\n")
    
    if (is.matrix(cov)) {
      SIGMA <- cov
      COR <- cov2cor(SIGMA)
    }
    else if (cov == "diag") {
      SIGMA <- diag(sdVAL^2, nrow = length(VARS), ncol = length(VARS))    
      colnames(SIGMA) <- rownames(SIGMA) <- colnames(DATA)
      COR <- cov2cor(SIGMA)
    }
    else if (cov == "full") {
      SIGMA <- cov(DATA, use = "pairwise.complete.obs")
      COR <- cor(DATA, use = "pairwise.complete.obs")
    }
  }  
  
  ## version 1.0-5: replace possible NA's/NaN's in covariance/correlation matrix with 0's
  SIGMA[is.na(SIGMA)] <- 0; SIGMA[is.nan(SIGMA)] <- 0
  COR[is.na(COR)] <- 0; COR[is.nan(COR)] <- 0
  
  ## version 1.0-5: No diagonals with 0, 
  if (any(diag(SIGMA) == 0)) {
    DIAG <- diag(SIGMA)
    DIAG[DIAG == 0] <- 1E-6
    diag(SIGMA) <- DIAG
  }
  
  ## This will bring the variables in 'data' and 'expr' in the same 
  ## order as in the covariance matrix
  m1 <- match(colnames(SIGMA), colnames(DATA))
  meanVAL <- meanVAL[m1]
  m2 <- match(colnames(SIGMA), VARS)
  
  ############ first- and second-order Taylor expansion-based error propagation ################
  ## first-order mean: eval(EXPR)
  ## version 1.0-4: continue with NA's when differentiation not possible
  MEAN1 <- try(eval(EXPR, envir = as.list(meanVAL)), silent = TRUE)
  if (!is.numeric(MEAN1)) {
    message("propagate: there was an error in calculating the first-order mean")
    MEAN1 <- NA
  }  
  
  ## evaluate gradient vector
  GRAD <- try(makeGrad(EXPR, m2), silent = TRUE)  
  if (!inherits(GRAD, "try-error")) evalGRAD <- try(sapply(GRAD, eval, envir = as.list(meanVAL)), silent = TRUE)
  if (inherits(GRAD, "try-error")) evalGRAD <- try(numGrad(EXPR, as.list(meanVAL)), silent = TRUE)  
  if (!inherits(evalGRAD, "try-error")) evalGRAD <- as.vector(evalGRAD) else evalGRAD <- NA
  
  ## first-order variance: g.S.t(g) 
  VAR1 <- try(as.numeric(t(evalGRAD) %*% SIGMA %*% matrix(evalGRAD)), silent = TRUE)  
  if (inherits(VAR1, "try-error")) {
    message("propagate: there was an error in calculating the first-order variance")
    VAR1 <- NA
  }
  
  ## second-order mean: firstMEAN + 0.5 * tr(H.S) 
  HESS <- try(makeHess(EXPR, m2), silent = TRUE)
  if (!inherits(HESS, "try-error")) evalHESS <- try(sapply(HESS, eval, envir = as.list(meanVAL)), silent = TRUE)
  if (inherits(HESS, "try-error")) evalHESS <- try(numHess(EXPR, as.list(meanVAL)), silent = TRUE)
  if (!inherits(evalHESS, "try-error")) evalHESS <- matrix(evalHESS, ncol = length(meanVAL), byrow = TRUE) else evalHESS <- NA  
 
  valMEAN2 <- try(0.5 * tr(evalHESS %*% SIGMA), silent = TRUE)
  if (!inherits(valMEAN2, "try-error")) {
    MEAN2 <- MEAN1 + valMEAN2
  } else {
    warning("propagate: there was an error in calculating the second-order mean")
    MEAN2 <- NA
  }
  
  ## second-order variance: firstVAR + 0.5 * tr(H.S.H.S)
  valVAR2 <- try(0.5 * tr(evalHESS %*% SIGMA %*% evalHESS %*% SIGMA), silent = TRUE)
  if (!inherits(valVAR2, "try-error")) {
    VAR2 <- VAR1 + valVAR2
  } else {
    message("propagate: there was an error in calculating the second-order variance")
    VAR2 <- NA
  }
  
  ## total mean and variance  
  if (exp.uc == "1st") {
    totalMEAN <- MEAN1
    totalVAR <- VAR1
  }
  if (exp.uc == "2nd") {
    totalMEAN <- MEAN2
    totalVAR <- VAR2 
  }
  errorPROP <- sqrt(totalVAR)  
  
  ## sensitivity index/contribution/relative contribution
  if (is.numeric(evalGRAD)) {
    sensitivity <- evalGRAD
    names(sensitivity) <- colnames(DATA)
    contribution <- outer(sensitivity, sensitivity, "*") * SIGMA
    rel.contribution <- abs(contribution)/sum(abs(contribution), na.rm = TRUE)
  } else sensitivity <- contribution <- rel.contribution <- NULL
  
  ## WS degrees of freedom, coverage factor and expanded uncertainty
  if (!is.null(dfVAL)) dfVAL[is.na(dfVAL)] <- 1000
  if (is.null(dfVAL)) dfVAL <- rep(1000, ncol(DATA))
  ws <- WelchSatter(ui = sqrt(diag(SIGMA)), ci = sensitivity, df = dfVAL, df.tot = df.tot, uc = errorPROP, alpha = alpha, k = k)
  
  ## confidence interval based on either first- or second-order mean
  if (exp.uc == "1st") confMEAN <- MEAN1 else confMEAN <- MEAN2
  confPROP <- confMEAN + c(-1, 1) * ws$u.exp
  names(confPROP) <- paste(c(alpha/2, 1 - alpha/2) * 100, "%", sep = "")
  
  ################# version 1.0-8 : Monte-Carlo simulation copula of t-distributions #####################
  ## if 'sim' data, don't create Monte Carlo data
  if (type == "raw" | type == "stat") {
    COPULA <- normalCopula(param = P2p(COR), dim = ncol(DATA), dispstr = "un")
    MARGINS <- rep("t", ncol(DATA))
    if (!is.null(df.tot)) dfVAL <- rep(df.tot, ncol(DATA))
    paramMargins <- lapply(dfVAL, function(df) list(df = df))
    MVD <- mvdc(copula = COPULA, margins = MARGINS, paramMargins = paramMargins)
    datSIM <- rMvdc(nsim, MVD); datSIM2 <- rMvdc(nsim, MVD)
    colnames(datSIM) <- colnames(datSIM2) <- colnames(DATA)
    varcor <- sqrt(dfVAL/(dfVAL - 2))
    datSIM <- sweep(datSIM, 2, varcor, "/"); datSIM2 <- sweep(datSIM2, 2, varcor, "/")
    datSIM <- sweep(datSIM, 2, sdVAL, "*"); datSIM2 <- sweep(datSIM2, 2, sdVAL, "*")
    datSIM <- sweep(datSIM, 2, meanVAL, "+"); datSIM2 <- sweep(datSIM2, 2, meanVAL, "+")
  } else {datSIM <- DATA; datSIM2 <- NULL}
  
  ## version 1.0-8: covariance and correlation of MC Copula matrix
  covMC <- cov(datSIM, use = "pairwise.complete.obs")
  corMC <- cor(datSIM, use = "pairwise.complete.obs")
  frobCOV <- norm(covMC - SIGMA, "F")/norm(SIGMA, "F")
  frobCOR <- norm(corMC - COR, "F")/norm(COR, "F")
  
  ## version 1.0-8: fit scaled/shifted t-distribution on Copula marginals
  if (check) {
    cat("Checking Copula margins with scaled/shifted t-distribution...\n")
    fitDIST <- apply(datSIM, 2, function(x) fitDistr(x, distsel = 5, nbin = 1000, verbose = FALSE, plot = "none")$bestpar)
    fitDIST[2, ] <- fitDIST[2, ] * sqrt(dfVAL/(dfVAL - 2))  # convert scale to sd
  } else fitDIST <- NULL
  
  ## try vectorized evaluation, which is much faster  
  resSIM <- try(eval(EXPR, envir = as.data.frame(datSIM)), silent = TRUE) 
  
  ## use 'row-wise' method if 'vectorized' throws an error
  if (inherits(resSIM, "try-error")) {
    message("propagate: using 'vectorized' evaluation gave an error. Switching to 'row-wise' evaluation...")
    resSIM <- apply(datSIM, 1, function(x) eval(EXPR, envir = as.list(x)))     
  }
    
  ## alpha-based confidence interval of MC simulations
  confSIM <- quantile(resSIM, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE) 
    
  ## warning in case of single evaluated result
  if (length(unique(resSIM)) == 1) message("propagate: Monte Carlo simulation gave unique repetitive values! Are all derivatives constants?")   
  
  outPROP <- c(mean.1 = MEAN1, mean.2 = MEAN2, u.1 = sqrt(VAR1), u.2 = sqrt(VAR2), 
               confPROP[1], confPROP[2])   
  
  outSIM <- c(mean.MC = mean(resSIM, na.rm = TRUE), u.MC = sd(resSIM, na.rm = TRUE), 
              median.MC = median(resSIM, na.rm = TRUE), mad.MC = mad(resSIM, na.rm = TRUE),
              confSIM[1], confSIM[2])
  
  outDAT <- rbind(meanVAL, sdVAL, dfVAL)
  
  OUT <- list(gradient = GRAD, evalGrad = evalGRAD,
              hessian = HESS, evalHess = evalHESS, 
              rel.contr  = rel.contribution, covMat = SIGMA, covMat.MC = covMC,
              corMat = COR, corMat.MC = corMC, frobCOV = frobCOV, frobCOR = frobCOR, 
              ws.df = floor(ws$ws.df), k = ws$k, u.exp = ws$u.exp, resSIM = resSIM, datSIM = datSIM, datSIM2 = datSIM2, 
              prop = outPROP, sim = outSIM, expr = EXPR, data = outDAT, alpha = alpha, nsim = nsim, check = fitDIST, type = type)
  
  class(OUT) <- "propagate"
  return(OUT)                                     
}

