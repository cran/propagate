colVars <- function(x, na.rm = TRUE)
{
  x <- as.matrix(x)
  NR <- nrow(x)
  NC <- ncol(x)
  N <- NR - colSums(is.na(x), na.rm = na.rm)
  CV <- (colSums(x^2, na.rm = na.rm) - N * colMeans(x, na.rm = na.rm)^2)/(N - 1) 
  CV[N < 2] <- NA
  CV
}

rowVars <- function(x, na.rm = TRUE)
{
  x <- as.matrix(x)
  NR <- nrow(x)
  NC <- ncol(x)
  N <- NC - rowSums(is.na(x), na.rm = na.rm)
  RV <- (rowSums(x^2, na.rm = na.rm) - N * rowMeans(x, na.rm = na.rm)^2)/(N - 1) 
  RV[N < 2] <- NA
  RV
}

colSDs <- function(x, na.rm = TRUE) sqrt(colVars(x, na.rm = na.rm))
rowSDs <- function(x, na.rm = TRUE) sqrt(rowVars(x, na.rm = na.rm))
