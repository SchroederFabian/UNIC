#' Expected Loss of the Threshold Classifier.
#' 
#' @description The function offers a method to select variables by univariate 
#' filtering based on the estimated loss of the optimal univariate threshold 
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
#' \item{etc}{ a numerical vector containing the etc values for every variable 
#' of dat.}
#' \item{p.val}{ the corresponding p-values of etc. (optional)}
#' \item{p.val.adj}{ the corresponding adjusted p-values of etc. (optional)}
#' 
#' @examples 
#' oc <- c(1, 3, 0.5)
#' class <- factor(c(rep(0, 25), rep(1, 25)), labels = c("neg", "pos"))
#' data <- data.frame("var1" = c(rnorm(25, 0, 1/2), rnorm(25, 1, 2)))
#' res <- etc(class, data, positive = "pos", oc, p.val = TRUE)
#' 
#' @export

etc <- function(class, data, positive = levels(class)[1], oc = c(1, 1, 0.5),  
                p.val = TRUE, plot = FALSE, adj.method = "BH") {
  
  # Error handling
  if (is.null(dim(data))) data <- as.data.frame(data)
  error_handling("etc", class, data, positive, oc, p.val, plot, adj.method)
  
  # Define variables
  pos <- class == positive
  neg <- !pos
  p <- ifelse(class(data) == "numeric", 1, ncol(data))
  if (class(data) == "numeric") data <- as.data.frame(data)
  npos <- as.numeric(table(class)[positive])
  nneg <- as.numeric(table(class)[levels(class)[!levels(class) == positive]])
  
  # Return opjects
  epe <- vector(mode = "numeric", length = p)
  names(epe) <- colnames(data)
  p.value <- vector(length = p)	
  
  # Calculate the prediction error for every variable in the data set
  for (i in 1:p) {
    
    ord <- order(data[, i])
    class.ord <- class[ord]
    feat.ord <- data[ord, i]
    
    fp1 <- c(0, cumsum(as.numeric(class.ord == positive)))
    fn1 <- c(nneg, nneg - cumsum(as.numeric(!class.ord == positive)))
    fp2 <- c(npos, npos - cumsum(as.numeric(class.ord == positive)))
    fn2 <- c(0, cumsum(as.numeric(!class.ord == positive)))
    
    epe[i] <- min(min(fp1 / npos * oc[2] * (1 - oc[3]) + fn1 / nneg * oc[1] * oc[3]), 
                  min(fp2 / npos * oc[2] * (1 - oc[3]) + fn2 / nneg * oc[1] * oc[3]))
  }
  
  # generate the null distribution and calculate the p-values
  if (p.val) {
    ND <- etc.genND(nneg, npos, oc[1], oc[2], oc[3])			
    p.val <- cumsum(as.numeric(ND$fav.perm / ND$pos.perm))[match(round(epe, 4), round(ND$val, 4))]
    p.value.adj <- p.adjust(p.val, method = adj.method)
    
    # Generate plot
    if (plot) {
      plot(stats::stepfun(ND$val, 
                          c(cumsum(as.numeric(ND$fav.perm / ND$pos.perm)), 1)), 
                          main = "Cumulative Null Distribution of ETC", 
                          xlab = "EPE", 
                          ylab = "", 
                          pch = ".")
      hist.info <- hist(epe, breaks = ND$val, plot = FALSE)
      points(hist.info$mids[hist.info$count != 0], 
             hist.info$counts[hist.info$count != 0] / sum(hist.info$counts[hist.info$count != 0]), 
             col = "red", 
             pch = 20)
      text(hist.info$mids[hist.info$count != 0], 
           y = hist.info$counts[hist.info$count != 0] / sum(hist.info$counts[hist.info$count != 0]), 
           labels = hist.info$counts[hist.info$count != 0], 
           col = "red")
    }
    
  } else { 
    p.value <- NULL
    p.value.adj <- NULL
    }
  
  return(list("etc" = epe, "p.value" = p.val, "p.value.adj" = p.value.adj))
  
}