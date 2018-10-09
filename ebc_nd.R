

ebc.genND <- function(class, data, positive, oc, equalVars = FALSE, robust = FALSE) {
  
  n1 <- sum(class == positive)
  n0 <- sum(class != positive)
 
  reps <- 100000
  data.h0 <- matrix(ncol = reps, nrow = length(class), data=rnorm(reps*length(class)))
  distr.h0 <- ebc(class = class, data = data.h0, oc = oc, positive = positive, p.val = FALSE, equalVars = equalVars, robust = robust)$ebc
  min.val <- min(distr.h0)
  
  # fit empirical cumulative distribution function
  h0.ecdf <- ecdf(distr.h0)
  
  # fit monotonic spline
  distr.h0.tail <- distr.h0[which(rank(distr.h0)<=200)]
  grid <- seq(0,max(distr.h0.tail), 0.01)
  pts <- which(h0.ecdf(grid)>0)
  h0.spline <- splinefun(c(0,grid[pts]), c(0,h0.ecdf(grid[pts])), method = "hyman")
  pval <- function(x) ifelse(x < min.val, h0.spline(x), h0.ecdf(x))

  return( pval )
  
}