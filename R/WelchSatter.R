WelchSatter <- function(ui, ci = NULL, df = NULL, df.tot = NULL, uc = NULL, alpha = 0.05, k = NULL)
{
  if (is.null(ci)) ci <- rep(1, length(ui))
  if (is.null(df)) stop("WelchSatter: Please supply 'df's!")
  if (length(ui) != length(df)) stop("WelchSatter: Different number of values in 'ui' and 'df'!")
  if (is.null(uc)) uc <- sqrt(sum((ci * ui)^2))
  ws.df <- (uc^4)/sum(((ci * ui)^4)/df)
  if (is.numeric(df.tot)) ws.df <- df.tot
  if (is.na(ws.df) | ws.df < 3) {
    cat("WARNING: Effective DOF < 3 indicates heavy-tailed distribution. Setting DOF = 3, but beware!\n")
    ws.df <- 3
  }
  if (is.null(k)) k <- qt(1 - alpha/2, floor(ws.df))
  u.exp <- k * uc
  return(list(ws.df = ws.df, k = k, u.exp = u.exp))   
}