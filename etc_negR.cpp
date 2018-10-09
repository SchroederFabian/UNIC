/* These functions calculate the number of favorable permutations of the null distribution
 * of the etc classifier for the negative right side for a given set of parameters fp, fn, 
 * n0, n1, c0, c1, pi0. negRght calculates the starting values and calls the recursive 
 * function negRght_rec. 
 */

#include <wrap.hpp>
#include <map>
#include <gmpxx.h>
#include <gmp.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

//[[Rcpp::plugins(cpp11)]]
void negRght_rec(std::map <std::vector<int>, mpz_class>& memo, mpz_class* sums, int level, int start, int stop, int& tp, int& fn, int& tn, int& fp, double& wght) {
  
  // memoization
  int arr[3] = {level, start, stop};
  std::vector<int> itm(arr, arr+3);
  
  if (memo.find(itm) != memo.end()) { // map lookup
    *sums = *sums + memo.at(itm);
  } else { // calculate key value
    if (level == 1) { // recursion base case
      *sums = *sums + start - stop + 1;
      memo.insert ( std::make_pair(itm, start - stop + 1) );
    } else {
      int level_new = level - 1;
      int mnpsr, mxpsr;
      double crt;
      crt = round((fn - level_new + 1) * wght *1e4)/1e4;
      mnpsr=ceil(crt);
      crt = round( (tn - fp - mnpsr - (level_new - tp) * wght) * 1e4 ) / 1e4;
      mxpsr = floor(crt);
      if (mxpsr >= 0) {
        int start_new = (fn + tn) - mnpsr - (fn - level_new);  
        int stop_new = std::max(level_new, tn - mnpsr + level_new - mxpsr);
        mpz_class bfr(*sums);
        for (int i = start; i >= stop; i--) {
          start_new = std::min(i-1, start_new);
          negRght_rec(memo, sums, level_new, start_new, stop_new, tp, fn, tn, fp, wght);
        } 
        mpz_class aftr(*sums);
        memo.insert ( std::make_pair(itm, aftr - bfr) );
        
      } 
    }
  }
}


// [[Rcpp::export]]
SEXP negRght(Rcpp::NumericMatrix clst, int n0, int n1, double c0, double c1, double pi0) {
  
  int ncol = clst.ncol();
  std::vector <mpz_class> perm(ncol);
  
  for (int i = 0; i < ncol; i = i + 1) { 
    
    int fp = clst(0,i);
    int fn = clst(1,i);
    int tp = n1 - fn;
    int tn = n0 - fp;
    int n = n0 + n1;
    mpz_class sums(0);
    double crt;
    
    std::map <std::vector<int>, mpz_class> memo; // create map
    
    if (fn == 0) {
      
      sums = 1;
      perm[i] = sums;
      
    } else {
      
      int level = fn;
      double wght = c1 / c0 * n0 / n1 * (1 - pi0) / pi0;
      crt = wght;
      int mnpsr = (fabs(ceil(crt) - crt) < 1e-9) ? round(crt) : ceil(crt);
      crt = (tn - fp - mnpsr - (level - tp) * wght);
      int mxpsr = (fabs(floor(crt) - crt) < 1e-9) ? round(crt) : floor(crt);

      if (mxpsr >= 0) {
        int start = (fn + tn) - mnpsr;  
        int stop = std::max(level, tn - mnpsr + level - mxpsr);
        
        // start recursion
        negRght_rec(memo, &sums, level, start, stop, tp, fn, tn, fp, wght);
        perm[i] = sums;
      }
    }
    memo.clear(); // delete map
  }
  return(wrap(perm));
}
