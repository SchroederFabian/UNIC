#' Generate the Null Distribution of EIC.
#' 
#' @description The function offers a heuristic aproach to generating the null 
#' distribution of the eic classifier under the null hypothesis that the 
#' distributions of the positive and the negative class are identical.
#' 
#' @param n0 an integer indicating the number of negative instances in the 
#' sample.
#' @param n1 an integer indicating the number of positive instances in the 
#' sample.
#' @param c0 the cost of misclassifying a negative instance.
#' @param c1 the cost of misclassifying a positive instance.
#' @param pi0 a real number between 0 and 1 indicating the percentage of 
#' negative instances in the population.
#' 
#' @return a list containing three components:
#' \item{val}{ a vector with the number of possible values than eic can take.}
#' \item{pos.perm}{ the number of possible permutations for the given number of 
#' positives and negatives.}
#' \item{fav.perm}{ a vector with the number of favorable permutations for 
#' every value in val.}
#' 
#' @examples 
#' eic.genND(25, 27, 1, 3, 0.5)
#' 
#' @export

eic.genND <- function(n0, n1, c0, c1, pi0, data, class, positive, plot = FALSE) {
  
  # generate list with pairs (fp, fn), for which to calculate the favorable permutations
  mat <- round( outer( seq(0, n0) / n0 * c0 * pi0, seq(0, n1) / n1 * c1 * (1 - pi0), FUN = "+"), 4)
  val <- sort(as.vector(mat)[!duplicated(as.vector(mat))])
  val <- val[val <= min(c1 * (1 - pi0), c0 * pi0)]
  
  # generate list with pairs (fp, fn), for which to calculate the favorable permutations
  clst <- list()
  for (v in val) { clst[[as.character(v)]] <- t(which(mat == v, arr.ind = TRUE) - 1) }
  
  # select only those cost values where all pairs of (fp, fn) do not exceed floor (min(n0, n1) / 2)
  sel <- sapply(clst, function(y) all(apply(y, 2, function(x) all(x <= floor(min(n0,n1) / 2 ) ))))
  clst <- clst[sel]
  
  # unlist clst
  clst.vec <- clst[[1]]
  colnames(clst.vec)[1] <- paste0(clst[[1]], collapse=",")
  for (i in 2:length(clst)) { 
    clst.vec <- cbind(clst.vec, clst[[i]])
    colnames(clst.vec)[(ncol(clst.vec) - ncol(clst[[i]]) + 1) : ncol(clst.vec)] <- apply(clst[[i]], 2, function(x) paste0(x, collapse=","))
  }
  
  # number of possible permutations
  pos.perm <- gmp::chooseZ(n0 + n1, n0) 
  
  # number of favorable permutations for every value
  fav.perm <- rep(gmp::as.bigz(0), length(val)) 

  # calculate the exact number of favorable permutation for the first floor(min(n1,n0)/2) values
  posInt.lst <- posInt(clst.vec, n0, n1, c0, c1, pi0)
  negInt.lst <- negInt(clst.vec, n0, n1, c0, c1, pi0)
  posComp.lst <- posComp(clst.vec, n0, n1, c0, c1, pi0)
  negComp.lst <- negComp(clst.vec, n0, n1, c0, c1, pi0)
  fav.perm.vec <- posInt.lst * negComp.lst + posComp.lst * negInt.lst
  
  cntr <- 1
  
  for (i in 1:(floor(min(n0, n1) / 2))){
    erg <- gmp::as.bigz(0)

    for (j in 1:ncol(clst[[i]])) {
      
     erg <- gmp::add.bigz(erg, fav.perm.vec[cntr])
      cntr <- cntr + 1
    }
    fav.perm[i] <- erg
  }

  ind.count <- which(fav.perm != 0)
  
  
  ## ESTIMATE THE NUMBER OF FAVORABLE PERMUTATIONS BY MEANS OF RANDOM PERMUTATIONS
  # for all values that exhibit more than 5% of overall permutations
  
  reps <- 20000  # number of random permutations
  fav.perm.rnd <- rep(0, length(val))
  names(fav.perm.rnd) <- val
  dt <- matrix(nrow = nrow(data), ncol = reps) 
  for (l in 1:ncol(dt)) { dt[,l] <- sample(data[,sample(ncol(data),1)]) } # fill the matrix with randomly permuted vectors
  cvals <- eic(class, dt, c(c0, c1, pi0), positive = positive, p.val = FALSE, plot = FALSE, adj.method = NULL)
    
  for (i in 1:length(val)) { fav.perm.rnd[i] <- sum(round(cvals$eic, 4) == round(val[i], 4) ) }
    
  ind.perm <- fav.perm.rnd >= 0.01 * reps
  fav.perm.rnd <- round(fav.perm.rnd * as.numeric(pos.perm / reps))
  fav.perm.rnd[!ind.perm] <- 0
  fav.perm.rnd[ind.count] <- 0
  fav.perm[fav.perm.rnd != 0] <- fav.perm.rnd[fav.perm.rnd != 0]

  
  ## ESTIMATE THE NUMBER OF FAVORABLE PERMUTATIONS FOR ALL OTHER VALUES BY MEANS OF INTRAPOLATION OR EXTRAPOLATION
  
  # for those values between the maximal value of the exact counting schema and the minimal value of the permutation
  if (max(ind.count) < min(which(ind.perm)) ) {
    
    val.r.imp <- c(ind.count[(length(ind.count) - 2) : length(ind.count)], which(ind.perm)[c(1,2,3)])
    val.2.imp <- (max(ind.count) + 1):(min(which(ind.perm))-1)
    ptsSpline.mid <- spline(val[val.r.imp], fav.perm[val.r.imp], xout = val[val.2.imp], method = "hyman")
    
    # plot
    #ptsSpline.plt <- spline(val[val.r.imp], fav.perm[val.r.imp], n = 200, method = "hyman")
    #plot(ptsSpline.plt, type = "l")   
    #points(val[val.r.imp], fav.perm[val.r.imp])
    #points(val[val.2.imp], ptsSpline.mid$y, col="red")

    fav.perm[val.2.imp] <- ptsSpline.mid$y

  }
 
  # for those values at the right end
 
  if (sum(fav.perm == 0) > 1) {
    
    if (sum(fav.perm) < pos.perm) {
      val.ext <- fav.perm == 0
      # intrapolate the cumulative density
      fav.perm.cum <- cumsum(fav.perm)
      fav.perm.cum[val.ext] <- NA
      fav.perm.cum[length(fav.perm.cum)] <- pos.perm 
      val.r.imp <- c(which(ind.perm)[(sum(ind.perm)-2):sum(ind.perm)], length(val))
      val.2.imp <- (max(which(ind.perm)) + 1):(length(val) - 1)
      ptsSpline.rght <- spline(val[val.r.imp], fav.perm.cum[val.r.imp], xout = val[val.2.imp], method = "hyman")
      fav.perm.cum[val.2.imp] <- round(ptsSpline.rght$y)
      
      dfs <- as.numeric(fav.perm.cum[c(val.2.imp, length(val))] - fav.perm.cum[c(val.2.imp[1] -1) : max(val.2.imp)])
      fav.perm[ (max(which(ind.perm)) + 1) : length(val)] <- sapply(dfs, function(z) max(0,z))
      
    } else {
      # round bigq not implemented yet
      #fav.perm[-ind.count] <- round(fav.perm[-ind.count] * (pos.perm - sum(fav.perm[ind.count])) / sum(fav.perm[-ind.count]))
    }
  } else if (sum(fav.perm == 0) == 1) {
    fav.perm[fav.perm == 0] <- max(pos.perm - sum(fav.perm), 0)
  }
  
#  if (plot) {
#    plot(stepfun(val, c(0, as.numeric(fav.perm / pos.perm) )), main = "Null Distribution of EIC")
#    points(val[ind.count], as.numeric(fav.perm / pos.perm)[ind.count], col="red")
#    points(val[ind.perm], as.numeric(fav.perm / pos.perm)[ind.perm], col="green")
#    points(val[-c(ind.count, which(ind.perm))], as.numeric(fav.perm / pos.perm)[-c(ind.count, which(ind.perm))], col="blue")
#    legend("topleft", legend = c("exact", "permutation", "intrapolation"), col = c("red", "green", "blue"), pch = rep(1,3))
#  }
 
  
  return(list("val" = val, "pos.perm" = pos.perm, "fav.perm" = fav.perm))
}  


