bigcor <- function(
    x, 
    y = NULL,
    fun = c("cor", "cov"), 
    size = 2000, 
    verbose = TRUE,
    file = "result.h5",
    ...)
{
  fun <- match.arg(fun)
  if (fun == "cor") FUN <- cor else FUN <- cov
  if (fun == "cor") STR <- "Correlation" else STR <- "Covariance" 
  if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")
  
  NCOL <- ncol(x)
  if (!is.null(y)) YCOL <- NCOL(y)
  
  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST  
  NBLOCKS <- NCOL %/% size
  
  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)  
  if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  
  ## initiate time counter
  timeINIT <- proc.time()[3]
  
  ## create HDF5 file and preallocate dataset
  if (file.exists(file)) file.remove(file)
  h5file <- hdf5r::H5File$new(file, mode = "w")
  
  # Use larger chunks aligned with block size and disable compression
  if (is.null(y)) {
    resMAT <- h5file$create_dataset(
      name = "matrix",
      dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
      space = hdf5r::H5S$new("simple", dims = c(NCOL, NCOL), maxdims = c(NCOL, NCOL)),
      chunk_dims = c(size, size)
    )
  } else {
    resMAT <- h5file$create_dataset(
      name = "matrix",
      dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
      space = hdf5r::H5S$new("simple", dims = c(NCOL, YCOL), maxdims = c(NCOL, YCOL)),
      chunk_dims = c(size, YCOL)
    )
  }
  
  ## iterate through each block combination
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]    
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]    
    
    timeNOW <- proc.time()[3] - timeINIT
    
    ## if y = NULL
    if (is.null(y)) {
      if (verbose) cat(sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... %.2f sec\n", i, STR, COMB[1], COMB[2], length(G1), length(G2), signif(timeNOW, 2)))
      RES <- FUN(x[, G1], x[, G2], ...)
      
      # Write operations
      resMAT[G1, G2] <- RES
      
      if (COMB[1] != COMB[2]) resMAT[G2, G1] <- t(RES)
    } else ## if y = smaller matrix or vector  
    {
      if (verbose) cat(sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... %.2f sec\n", i, STR, COMB[1], length(G1), YCOL, signif(timeNOW, 2)))    
      RES <- FUN(x[, G1], y, ...)
      resMAT[G1, ] <- RES             
    }
    
    gc()
  } 
  
  # Close the HDF5 file
  h5file$close_all()
  if (verbose) cat("=> Results saved to '", file, "'")
  
  return(invisible(file))
}
