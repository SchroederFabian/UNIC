/* These functions calculate the number of favorable permutations of the null distribution
 * of the eic classifier for a negative interval for a given set of parameters fp, fn, 
 * n0, n1, c0, c1, pi0. negInt calculates the starting values and calls the recursive 
 * function negInt_rec. 
 */

#include <wrap.hpp>
#include <map>
#include <gmpxx.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

void negInt_rec(std::map <std::vector<int>, mpz_class>& memo_negInt, mpz_class* sum_negInt, int level, int start, int stop, int fp, int fn, int n0, int n1, double c0, double c1, double pi0) {
  
  int arr[5] = {level, start, stop, fp, fn}; // map key
  std::vector<int> itm(arr, arr+5);
  
  if (memo_negInt.find(itm) != memo_negInt.end()) { // map look up
    
    *sum_negInt = *sum_negInt + memo_negInt.at(itm);
    
  } else { // if key does not exist in map, calculate the corresponding value
    
    if (level == 1) { // end of recursion
      
      *sum_negInt = *sum_negInt + stop - start + 1;
      
    } else { // recursion step
      
      // calculate new recursion arguments
      int level_new = level - 1;
      double crt = (fn - level_new + 1) * c1 / c0 * n0 / n1;
      int mnpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
      crt = level_new * c1 / c0 * n0 / n1;
      int mxpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
      int start_new = mnpsr + (fn - level_new) + 1;
      int stop_new = (n0 - fp + fn) - (mxpsr + (level_new - 1));
      
      mpz_class bfr(*sum_negInt);
      
      for (int i = start; i <= stop; i = i + 1) {
        // call recursion
        negInt_rec(memo_negInt, sum_negInt, level_new, std::max(i+1, start_new), stop_new, fp, fn, n0, n1, c0, c1, pi0 );   
      }
      
      mpz_class aftr(*sum_negInt);
      memo_negInt.insert( std::make_pair(itm, aftr - bfr)); // save pair of (key, value) into map
      
    }
  }
}


// [[Rcpp::export]]
SEXP negInt(Rcpp::NumericMatrix clst, int n0, int n1, double c0, double c1, double pi0) {
  
  int ncol = clst.ncol();
  std::vector <mpz_class> perm(ncol);
  std::map <std::vector<int>, mpz_class> memo; // create map
  
  for (int i = 0; i < ncol; i = i + 1) {
    
    mpz_class sums(0);
    int fp = clst(0,i);
    int fn = clst(1,i);
    
    if (fn == 0) {
      sums = 1;
      perm[i] = sums;
    } else {
      
      // calculate initial values of recursion
      int level_init = fn;
      double crt = c1 / c0 * n0 / n1;
      int mnpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
      crt = fn * c1 / c0 * n0 / n1;
      int mxpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
      
      int start_init = mnpsr + 1;
      int stop_init = (n0 - fp + fn) - (mxpsr + (fn - 1));
      
      if (start_init > stop_init) {
        sums = 0;
        perm[i] = sums;
      } else {
        
        // start recursion
        negInt_rec(memo, &sums, level_init, start_init, stop_init, clst(0,i), clst(1,i), n0, n1, c0, c1, pi0);
        perm[i] = sums;
      
      }
    }
  }
  
  memo.clear(); // delete map
  return(wrap(perm));
}
