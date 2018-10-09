#' Generate the Null Distribution of the Threshold Classifier.
#' 
#' @description The function offers an algorithmic aproach to generating the 
#' null distribution of the etc classifier under the null hypothesis that the 
#' distributions of the positive and the negative class are identical.
#' 
#' @param n0 integer indicating the number of negative instances in the 
#' sample.
#' @param n1 integer indicating the number of positive instances in the 
#' sample.
#' @param c0 the cost of misclassifying a negative instance.
#' @param c1 the cost of misclassifying a positive instance.
#' @param pi0 a real number between 0 and 1 indicating the percentage of 
#' negative instances in the population.
#' 
#' @return a list containing three components:
#' \item{val}{ a vector with the number of possible values than etc can take.}
#' \item{pos.perm}{ the number of possible permutations for the given number of 
#' positives and negatives.}
#' \item{fav.perm}{ a vector with the number of favorable permutations for every 
#' value in val.}
#' 
#' @examples 
#' etc.genND(25, 27, 1, 3, 0.5)
#' 
#' @export

etc.genND <- function(n0, n1, c0, c1, pi0) {
  
  # error handlig
  
  
  # generate list with pairs (fp, fn), for which to calculate the number of 
  # favorable permutations
  mat <- round(outer(seq(0, n0) / n0 * c0 * pi0, seq(0, n1) / n1 * c1 * (1 - pi0), FUN = "+"), 4)
  val <- sort(as.vector(mat)[!duplicated(as.vector(mat))])
  val <- val[val <= min(c1 * (1 - pi0), c0 * pi0)]
  clst <- list()
  for (v in val) { clst[[as.character(v)]] <- t(which(mat == v, arr.ind = TRUE) - 1) }
  
  # unlist clst
  clst.vec <- clst[[1]]
  colnames(clst.vec)[1] <- paste0(clst[[1]], collapse = ",")
  for (i in 2:length(clst)) { 
    clst.vec <- cbind(clst.vec, clst[[i]])
    colnames(clst.vec)[(ncol(clst.vec) - ncol(clst[[i]]) + 1) : ncol(clst.vec)] <- 
      apply(clst[[i]], 2, function(x) paste0(x, collapse = ","))
  }
  
  # calculate the number of possible permutations
  pos.perm <- gmp::chooseZ(n0 + n1, n0)

  # for every pair in clst calculate the number of favorable permutations
  fav.perm <- rep(gmp::as.bigz(0), length(val))
  
  posLeft.lst <- posLeft(clst.vec, n0, n1, c0, c1, pi0)
  negRght.lst <- negRght(clst.vec, n0, n1, c0, c1, pi0)
  posRght.lst <- posRght(clst.vec, n0, n1, c0, c1, pi0)
  negLeft.lst <- negLeft(clst.vec, n0, n1, c0, c1, pi0)
  
  fav.perm.vec <- posLeft.lst * negRght.lst + negLeft.lst * posRght.lst
  
  # add all permutations of pair (fp, fn) that add up to the same cost
  cntr <- 1
  for (i in 1:length(clst)){
    erg <- gmp::as.bigz(0)
    
    for (j in 1:ncol(clst[[i]])) {
      
      erg <- gmp::add.bigz(erg, fav.perm.vec[cntr])
      cntr <- cntr + 1
    }
    fav.perm[i] <- erg
  }
  
  return(list("val" = val, "pos.perm" = pos.perm, "fav.perm" = fav.perm))
}
