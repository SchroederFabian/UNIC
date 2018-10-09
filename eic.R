#' Expected Loss of the Interval Classifier.
#' 
#' @description The function offers a method to select variables by univariate
#' filtering based on the estimated loss of the optimal univariate interval 
#' classifer. No parametric assumption about the class conditional 
#' distributions is required.
#' 
#' @param class a factor vector indicating the class membership of the 
#' instances. Must have exactly two levels.
#' @param data a data frame with variables in columns.
#' @param oc a vector containing three elements. oc[1], the cost of 
#' misclassifying a negative instance, oc[2], the cost of missclassifying a
#' positive instance, and oc[3], the share of negative instances in the 
#' population.
#' @param positive a character object indicating the factor label of the 
#' positive class.
#' @param p.val a logical indicating whether the p-values of etc values under 
#' the null hypothesis that both classes are equal should be calculated. The 
#' exact null distribution is calculated by means of a recursive algorithm.
#' @param adj.method a character string indicating the method with which to 
#' correct the p-values for multiple testing. See ?p.adjust.
#' @param plot a logical. If TRUE a plot of the null distribution will be 
#' generated.
#' 
#' @return a list containing three components:
#' \item{eic}{ a numerical vector containing the etc values for every variable 
#' of dat.}
#' \item{p.val}{ the corresponding p-values of etc. (optional)}
#' \item{p.val.adj}{ the corresponding adjusted p-values of etc. (optional)}
#' 
#' @examples 
#' oc <- c(1,3,0.5)
#' class <- factor(c(rep(0, 25), rep(1, 25)), labels=c("neg", "pos"))
#' data <- data.frame("var1"=c(rnorm(25, 0, 1/2), rnorm(25, 1, 2)))
#' res <- eic(class, data, oc, positive="pos", p.val=TRUE)
#' 
#' @export


eic <- function(class, data, positive = levels(class)[1], oc = c(1, 1, 0.5),  
                p.val = TRUE, plot = FALSE, adj.method = "BH") {
  
  # Error handling
  if (is.null(dim(data))) data <- as.data.frame(data)
  error_handling("eic", class, data, positive, oc, p.val, plot, adj.method)
  
  # Define variables
  p <- ifelse(class(data) == "numeric", 1, ncol(data))
  n1 <- as.numeric(table(class)[positive])
  n0 <- length(class) - n1
  n <- n0 + n1
  d1 <- - oc[2] / n1 * (1 - oc[3])
  d0 <- oc[1] / n0 * oc[3]
  
  # Return objects
  epe <- vector(mode="numeric", length = p)
  names(epe) <- colnames(data)
  p.value <- vector(length = p)
  
  # Calculate the prediction error for every variable in the data set
  
  # order the class labels according to the ranks of the variables
  data.ord <- matrix(nrow = length(class), data = as.numeric(class)[apply(data, 2, order)])
  
  for (i in 1:p) {
    min.val <- min(oc[2] * (1 - oc[3]), oc[1] * oc[3])
    
    for (j in 1:n) { # positive interval
      row.min <- min(cumsum( c(oc[2] * (1 - oc[3]), c(d0,d1)[data.ord[j:n, i]]) ))
      min.val <- min(min.val, row.min)
    }
    
    for (j in 1:n) { # negative interval
      row.min <- min(cumsum(c(oc[1] * oc[3], c(-d0, -d1)[data.ord[j:n, i]] ) ))
      min.val <- min(min.val, row.min)
    }
    
    epe[i] <- min.val
    
  }
  
  # Generate the null distribution and calculate the p-values
  
  if (p.val) {
    ND <- eic.genND(n0, n1, oc[1], oc[2], oc[3], data, class, positive)			
    p.val <- cumsum(as.numeric(ND$fav.perm/ND$pos.perm))[match(round(epe,4), round(ND$val,4))]
    
    if (plot) {
      plot(stats::stepfun(ND$val, c(cumsum(as.numeric(ND$fav.perm / ND$pos.perm)), 1)), 
           main = "Cumulative Null Distribution of ETC", xlab = "EPE", ylab = "", pch = ".")
      hist.info <- hist(epe, breaks = ND$val, plot = FALSE)
      points(hist.info$mids[hist.info$count != 0], hist.info$counts[hist.info$count != 0] / sum(hist.info$counts[hist.info$count != 0]), col = "red", pch = 20)
      text(hist.info$mids[hist.info$count != 0], y = hist.info$counts[hist.info$count != 0] / sum(hist.info$counts[hist.info$count != 0]), labels = hist.info$counts[hist.info$count != 0], col = "red")
      
      p.val.adj <- ifelse(is.null(adj.method), NULL, p.adjust(p.val, method = adj.method))
    }
    
  } else { 
    p.val <- NULL 
    p.val.adj <- NULL
  }
  
  return(list("eic" = epe, "p.value" = p.val, "p.value.adj" = p.val.adj))

}



  
  
  
