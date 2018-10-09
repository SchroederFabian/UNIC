/* These functions calculate the number of favorable permutations of the null distribution
 * of the etc classifier for the positive left side for a given set of parameters fp, fn, 
 * n0, n1, c0, c1, pi0. posLeft calculates the starting values and calls the recursive 
 * function posLeft_rec. 
 */

#include <wrap.hpp>
#include <map>
#include <gmpxx.h>
#include <gmp.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

//[[Rcpp::plugins(cpp11)]]
void posLeft_rec(std::map <std::vector<int>, mpz_class>& memo, mpz_class* sums, int level, int start, int stop, int& tp, int& fn, int& tn, int& fp, double& wght) {
  
  // memoization
  int arr[3] = {level, start, stop};
  std::vector<int> itm(arr, arr+3);
  
  if (memo.find(itm) != memo.end()) { // map lookup
    *sums = *sums + memo.at(itm);
  } else {
    if (level == 1) { // recursion base case
      *sums = *sums + start - stop + 1;
      memo.insert ( std::make_pair(itm, start - stop + 1) );
    } else {
      
      int level_new = level - 1;
      int mnpsr, mxpsr;
      double crt;
      crt = round((fp - level_new + 1) * wght *1e4)/1e4;
      mnpsr = (fabs(ceil(crt) - crt) < 1e-6) ? round(crt + 1) : ceil(crt);
      crt = round((tp - fn - mnpsr - (level_new - tn) * wght) * 1e4) / 1e4;
      mxpsr = floor(crt);
      if (mxpsr >= 0) {
        int start_new = (fp + tp) - mnpsr - (fp - level_new);  
        int stop_new = std::max(level_new, tp - mnpsr + level_new - mxpsr);
        mpz_class bfr(*sums);
        for (int i = start; i >= stop; i--) {
          start_new = std::min(i-1, start_new);
          posLeft_rec(memo, sums, level_new, start_new, stop_new, tp, fn, tn, fp, wght);
        }
        mpz_class aftr(*sums);
        memo.insert ( std::make_pair(itm, aftr - bfr) );
      } 
    }
  }
}


// [[Rcpp::export]]
SEXP posLeft(Rcpp::NumericMatrix clst, int n0, int n1, double c0, double c1, double pi0) {
  
  int ncol = clst.ncol();
  std::vector <mpz_class> perm(ncol);
  
  for (int i = 0; i < ncol; i = i + 1) { 
  
  int fp = clst(0,i);
  int fn = clst(1,i);
  int tp = n1 - fn;
  int tn = n0 - fp;
  int n = n0 + n1;
  mpz_class sums(0);
  
  std::map <std::vector<int>, mpz_class> memo; // create map
  
  if (fp == 0) {
    
    sums = (tp == 0) ? 0 : 1;
    perm[i] = sums;
    
  } else {
    
    int level = fp;
    double wght = c0 / c1 * n1 / n0 * pi0 / (1 - pi0);
    double crt = wght;
    int mnpsr = (fabs(ceil(crt) - crt) < 1e-6) ? round(crt+1) : ceil(crt);
    crt = round((tp - fn - mnpsr - (level - tn)*wght)*1e4)/1e4;
    int mxpsr = floor(crt);

    if (mxpsr >= 0) {
      int start = fp + tp - mnpsr;
      int stop = std::max(level, tp - mnpsr + level - mxpsr);
      
      // start recursion
      posLeft_rec(memo, &sums, level, start, stop, tp, fn, tn, fp, wght);
      perm[i] = sums;
    }
  }
  memo.clear(); // delete map
  }
  return(wrap(perm));
}
  