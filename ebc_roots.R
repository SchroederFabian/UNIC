
roots <- function(mu0, mu1, sigma0, sigma1, c0, c1, pi0) {
  
  a <- 1 / (2 * sigma0^2) - 1 / (2 * sigma1^2)
  b <- mu1 / sigma1^2 - mu0 / sigma0^2
  c <- mu0^2 / (2 * sigma0^2) - mu1^2 / (2 * sigma1^2) + log(sigma0 / sigma1) - log(pi0 / (1 - pi0) * c0 / c1)
  
  if (b^2 - 4 * a * c >= 0) { return(c((- b - sqrt(b^2 - 4 * a * c)) / (2 * a), (- b + sqrt(b^2 - 4 * a * c)) / (2 * a)))										
  } else {
    return(c(NA, NA))
    }
}
