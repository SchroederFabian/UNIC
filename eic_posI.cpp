/* These functions calculate the number of favorable permutations of the null distribution
 * of the eic classifier for a positive interval for a given set of parameters fp, fn, 
 * n0, n1, c0, c1, pi0. posInt calculates the starting values and calls the recursive 
 * function posInt_rec. 
 */

#include <wrap.hpp>
#include <map>
#include <gmpxx.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

void posInt_rec(std::map <std::vector<int>, mpz_class>& memo, mpz_class* sum_posInt, int level, int start, int stop, int fp, int fn, int n0, int n1, double c0, double c1, double pi0) {
  
  int arr[5] = {level, start, stop, fp, fn}; // map key
  std::vector<int> itm(arr, arr+5);
  
  if (memo.find(itm) != memo.end()) { // map look up
    
    *sum_posInt = *sum_posInt + memo.at(itm);
    
  } else { // if key does not exist in map, calculate the corresponding value
    
    if (level == 1) { // end of recursion
      
      *sum_posInt = *sum_posInt + stop - start + 1;
      
    } else { // recursion step
      
      // calculate new recursion arguments
      int level_new = level - 1;
      double crt = (fp - level_new + 1) * c0 / c1 * n1 / n0;
      int mnpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
      crt = level_new * c0 / c1 * n1 / n0;
      int mxpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
      int start_new = mnpsr + (fp - level_new) + 1;
      int stop_new = (n1 - fn + fp) - (mxpsr + (level_new - 1));
      
      mpz_class bfr(*sum_posInt);
      
      for (int i = start; i <= stop; i = i + 1) {
        // call recursion
        posInt_rec(memo, sum_posInt, level_new, std::max(i+1, start_new), stop_new, fp, fn, n0, n1, c0, c1, pi0 );   
      }
      
      mpz_class aftr(*sum_posInt);
      memo.insert( std::make_pair(itm, aftr - bfr)); // save pair of (key, value) into map
      
    }
  }
}


// [[Rcpp::export]]
SEXP posInt(Rcpp::NumericMatrix clst, int n0, int n1, double c0, double c1, double pi0) {
  
  int ncol = clst.ncol();
  std::vector <mpz_class> perm(ncol);
  std::map <std::vector<int>, mpz_class> memo; // create map
  
  for (int i = 0; i < ncol; i = i + 1) {
    
    mpz_class sums(0);
    int fp = clst(0,i);
    int fn = clst(1,i);

    if (fp == 0) {
      sums = 1;
      perm[i] = sums;
    } else {
      
      // calculate initial values for recursion
      int level_init = fp;
      double crt = c0 / c1 * n1 / n0;
      int mnpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
      crt = fp * c0 / c1 * n1 / n0;
      int mxpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
      int start_init = mnpsr + 1;
      int stop_init = (n1 - fn + fp) - (mxpsr + (fp - 1));
      
      // start recursion
      posInt_rec(memo, &sums, level_init, start_init, stop_init, clst(0,i), clst(1,i), n0, n1, c0, c1, pi0);
      perm[i] = sums;
      
    }
  }
  
  memo.clear(); // delete map
  return(wrap(perm));
}

