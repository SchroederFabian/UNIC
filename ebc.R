#' Expected Loss of the Bayesian Classifier
#' 
#' @description The function offers a method to select variables by univariate
#' filtering based on the estimated loss of the univariate Bayesian Classifer.
#' The statistic requires the parametric assumption that the variable consists
#' of a mixture of Gaussian variables.
#' 
#' @param class a factor vector indicating the class membership of the 
#' instances. Must have exactly two levels.
#' @param data a data frame with the variables to filer in columns.
#' @param oc a vector containing three elements. oc[1], the cost of 
#' misclassifying a negative instance, oc[2], the cost of missclassifying a 
#' positive instance, and oc[3], the share of negative instances in the 
#' population.
#' @param positive a character object indicating the factor label of the 
#' positive class.
#' @param robust a logical indicating whether a robust estimator of the mean 
#' and variance of the two classes should be used.
#' @param p.val a logical indicating whether the p-values of ebc values under
#' the null hypothesis that both classes are equal should be calculated. 
#' Currently the null distribution is calculated by permutation.
#' @param adj.method a character string indicating the method with which to 
#' correct the p-values for multiple testing. See ?p.adjust.
#' 
#' @return a list containing three components:
#' \item{ebc}{ a numerical vector containing the etc values for every variable 
#' of dat.}
#' \item{p.val}{ the corresponding p-values of etc. (optional)}
#' \item{p.val.adj}{ the corresponding adjusted p-values of ebc. (optional)}
#' 
#' @examples 
#' oc <- c(1,3,0.5)
#' class <- factor(c(rep(0, 25), rep(1, 25)), labels=c("neg", "pos"))
#' data <- data.frame("var1"=c(rnorm(25, 0, 1/2), rnorm(25, 1, 2)))
#' res <- ebc(class, data, oc, positive="pos", p.val=TRUE)
#' @export

ebc <- function(class, data, positive = levels(class)[1], oc = c(1, 1, 0.5),  
                robust = FALSE,  p.val = TRUE, adj.method = "BH", equalVars = FALSE) {
  
  # Error handling
  if (is.null(dim(data))) data <- as.data.frame(data)
  error_handling("ebc", class, data, positive, oc, p.val, adj.method)
  
  # Define variables
  p <- ifelse(class(data) == "numeric", 1, ncol(data))
  n <- length(class)
  pos <- class == positive
  neg <- !pos
  
  # Return objects
  epe <- vector(mode = "numeric", length = p)
  names(epe) <- colnames(data)
  p.value <- vector(length = p)
  
  # Calculate the prediction error for every variable in the data set
  for (i in 1:p) {
    mu1 <- ifelse(robust, median(data[pos, i], na.rm = TRUE), 
                          mean(data[pos, i], na.rm = TRUE))
    mu0 <- ifelse(robust, median(data[neg, i], na.rm = TRUE), 
                          mean(data[neg, i], na.rm = TRUE))
    sigma1 <- ifelse(robust, mad(data[pos, i], na.rm = TRUE), 
                             sd(data[pos, i], na.rm = TRUE))
    sigma0 <- ifelse(robust, mad(data[neg, i], na.rm = TRUE), 
                             sd(data[neg, i], na.rm = TRUE))
    
    if (sigma1 == 0) sigma1 <- 1e-12
    if (sigma0 == 0) sigma0 <- 1e-12
    
    
    if (sigma0 == sigma1 | equalVars){ # equal variances
      
      d.std <- data[,i]
      d.std[pos] <- d.std[pos] - mean(data[pos,i], na.rm = TRUE)
      d.std[neg] <- d.std[neg] - mean(data[neg,i], na.rm = TRUE)
      sigma <- ifelse(robust, mad(d.std, na.rm = TRUE), sd(d.std, na.rm = TRUE))
      
      x <- log(oc[1] / oc[2] * oc[3] / (1 - oc[3])) * sigma^2 / (mu1 - mu0) + (mu1 + mu0)/2
      
      epe[i] <- ifelse(mu0 < mu1, 
                       oc[3] * oc[1] * (1 - pnorm(x, mu0, sigma)) + 
                         (1 - oc[3]) * oc[2] * pnorm(x, mu1, sigma), 
                       oc[3] * oc[1] * pnorm(x, mu0, sigma) + 
                         (1 - oc[3]) * oc[2] * (1 - pnorm(x, mu1, sigma)))
       
    } else { # unequal variances 
      
      if (mu0 < mu1) {
        xm <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[2]
        xa <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[1]
        if (is.na(xm)) { 
          epe[i] <- (min((1 - oc[3]) * oc[2], oc[3] * oc[1])) 
        } else {
          epe[i] <- ifelse (sigma0 < sigma1, 	
                            oc[3] * oc[1] * (1 - pnorm(xm, mu0, sigma0) + pnorm(xa, mu0, sigma0)) + 
                              (1 - oc[3]) * oc[2] * (pnorm(xm, mu1, sigma1) - pnorm(xa, mu1, sigma1)), 
                            oc[3] * oc[1] * (pnorm(xa, mu0, sigma0) - pnorm(xm, mu0, sigma0)) + 
                              (1 - oc[3]) * oc[2] * (1 - pnorm(xa, mu1, sigma1) + pnorm(xm, mu1, sigma1))) 
        }	
      } else {
        xm <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[1]
        xa <- roots(mu0, mu1, sigma0, sigma1, oc[1], oc[2], oc[3])[2]
        if (is.na(xm)) { 
          epe[i] <- (min((1 - oc[3]) * oc[2], oc[3] * oc[1])) 
        } else {
          epe[i] <- ifelse (sigma0 < sigma1, 	
                            oc[3] * oc[1] * (1 - pnorm(xa, mu0, sigma0) + pnorm(xm, mu0, sigma0)) + 
                            (1 - oc[3]) * oc[2] * (pnorm(xa, mu1, sigma1) - pnorm(xm, mu1, sigma1)), 
                            oc[3] * oc[1] * (pnorm(xm, mu0, sigma0) - pnorm(xa, mu0, sigma0)) + 
                            (1 - oc[3]) * oc[2] * (1 - pnorm(xm, mu1, sigma1) + pnorm(xa, mu1, sigma1)))
        }
      }
      if (sigma0 == 0 | sigma1 == 0) { 
        epe[i] <- NA 
        }
      if (is.na(xm)) {
        epe[i] <- (min((1 - oc[3]) * oc[2], oc[3] * oc[1])) 
        }
    }
  }
  

  # generate the null distribution and calculate the p-values
  if (p.val) {
   
    reps <- 100000
    data.h0 <- matrix(ncol = reps, nrow = length(class), data=rnorm(reps*length(class)))
    distr.h0 <- ebc(class = class, data = data.h0, oc = oc, positive = positive, p.val = FALSE)$ebc
    min.val <- min(distr.h0)
    
    # fit empirical cumulative distribution function
    h0.ecdf <- ecdf(distr.h0)
    
    # fit monotonic spline
    distr.h0.tail <- distr.h0[which(rank(distr.h0)<=200)]
    grid <- seq(0,max(distr.h0.tail), 0.01)
    pts <- which(h0.ecdf(grid)>0)
    h0.spline <- splinefun(c(0,grid[pts]), c(0,h0.ecdf(grid[pts])), method = "hyman")
    pval <- function(x) ifelse(x < min.val, h0.spline(x), h0.ecdf(x))
    p.value <- sapply(epe, pval)
    p.val.adj <- p.adjust(p.value, method = adj.method)
  } else { 
    p.value <- NULL
    p.val.adj <- NULL
    }
  
  return(list("ebc" = epe, "p.value" = p.value, "p.value.adj" = p.val.adj))
  
}
