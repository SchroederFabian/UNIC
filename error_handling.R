
error_handling <- function(func, class, data, positive, oc, p.val, plot, adj.method) {
  
  if (is.null(dim(data))) data <- as.data.frame(data)
  
  if (nlevels(factor(class)) != 2) stop("Class must consist of two groups")
  if (!positive %in% levels(class)) stop("Positive class doesn't exist")
  if (length(class) != dim(data)[1]) {
    stop("Dimensions of class and data do not correspond. Verify that features 
         are in columns.")
  }
  
  if (length(oc) != 3 | oc[3] < 0 | oc[3] > 1 | oc[1] < 0 | oc[2] < 0) {
    stop("The operating coditions do not meet the requirements.")
  }
  
  if (plot == TRUE & p.val == FALSE) {
    stop("In order to plot the result p.val must be set to TRUE.")
  }
  
  if (levels(class)[2] != positive) { 
    class <- factor(class, levels = levels(class)[c(2,1)] ) 
  }
  
  
}