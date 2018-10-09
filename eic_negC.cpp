/* These functions calculate the number of favorable permutations of the null distribution
 * of the eic classifier for a negative complement for a given set of parameters fp, fn, 
 * n0, n1, c0, c1, pi0. negComp calculates the starting values and calls the recursive 
 * function negComp_rec. 
 */
#include <wrap.hpp>
#include <map>
#include <gmpxx.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

void negComp_rec(std::map <std::vector<int>, mpz_class>& memo, mpz_class* sums, int level, int start, int stop, int tn_h, int tn_o, int fn_h, int fn_o, int n0, int n1, double c0, double c1, double pi0) {
  
  int arr[7] = {level, start, stop, tn_h, tn_o, fn_h, fn_o}; // map key
  std::vector<int> itm(arr, arr+7);
  
  if (memo.find(itm) != memo.end()) { // map look up
    
    *sums = *sums + memo.at(itm);
    
  } else { // if key does not exist in map, calculate the corresponding value
    
    if (level == 1) { // end of recursion
      
      *sums = *sums +  stop - start + 1;
      
    } else { // recursion step
      
      // calculate new recursion arguments
      int level_new = level - 1;
      double crt = level_new * c1 / c0 * n0 / n1;
      int mnpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
      crt = std::max((fn_o + (fn_h - level_new + 1)) * c1 / c0 * n0 / n1 - tn_o, -0.0001);
      int mxpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
      int start_new = mnpsr + level_new;
      int stop_new = (fn_h + tn_h) - (mxpsr + (fn_h - level_new));
      
      mpz_class bfr(*sums);
      
      for (int i = start; i <= stop; i = i + 1) {
        // call recursion
        negComp_rec(memo, sums, level_new, start_new, std::min(stop_new, i - 1), tn_h, tn_o, fn_h, fn_o, n0, n1, c0, c1, pi0);
      }
      
      mpz_class aftr(*sums);
      memo.insert( std::make_pair(itm, aftr - bfr)); // save pair of (key, value) into map
      
    }
  }
}


// [[Rcpp::export]]
SEXP negComp(Rcpp::NumericMatrix clst, int n0, int n1, double c0, double c1, double pi0) {
  
  // declaration of variables
  int ncol = clst.ncol();
  std::vector <mpz_class> perm(ncol); 
  std::map <std::vector<int>, mpz_class> memo; // create map
  int fp, fn, mid;
  int l_minpsr, r_minpsr, l_mxpsr, r_mxpsr;
  int l_level_init, r_level_init, l_start_init, r_start_init, l_stop_init, r_stop_init;
  double crt;
  mpz_class sum_lft(0), sum_rght(0), sum_tot(0);
  
  for (int i = 0; i < ncol; i = i + 1) {
    
    // setting
    fp = clst(0,i);
    fn = clst(1,i);
    sum_lft = 0;
    sum_rght = 0;
    sum_tot = 0;
    
    crt = static_cast<double>((n0 - fp + fn)/2);
    int mx = static_cast<int>(ceil(crt));
    mid = ((n0 - fp + fn) % 2 == 0) ? mx : -1;
    
    for (int l = 0; l <= mx; l = l + 1) {
      for (int fn_l = 0; fn_l <= fn; fn_l = fn_l + 1) {
        
        int tn_l = l - fn_l;  // true negatives on the left side
        int tn_r = (n0 - fp) - tn_l; // true negatives on the right side
        int fn_r = fn - fn_l; // false negatives on the right side
        
        sum_lft = 0;
        sum_rght = 0;
        
        // set l_minpsr and r_minpsr
        if (fn_l == 0) { l_minpsr = -1;
        } else {
          crt = fn_l * c1 / c0 * n0 / n1;
          l_minpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
        }
  
        if (fn_r == 0) { r_minpsr = -1;
        } else {
          crt = fn_r * c1 / c0 * n0 / n1;
          r_minpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
        }
  
        if (tn_l >= l_minpsr & tn_r >= r_minpsr) {
          
          // permutations of the left complement
          if (fn_l == 0) { sum_lft = 1;
          } else {
            
            // calculate initial recursion arguments
            l_level_init = fn_l;
            crt = (fn_r + 1) * c1 / c0 * n0 / n1 - tn_r;
            l_mxpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
            l_mxpsr = std::max(0, l_mxpsr);
            l_start_init = l_minpsr + (fn_l - 1) + 1;
            l_stop_init = (fn_l + tn_l) - l_mxpsr;
            
            if (l_start_init <= l_stop_init & l_stop_init <= fn_l + tn_l) {
              
              // start recursion
              negComp_rec(memo, &sum_lft, l_level_init, l_start_init, l_stop_init, tn_l, tn_r, fn_l, fn_r, n0, n1, c0, c1, pi0);

            } else { sum_lft = 0; }
          }
          
          // permutations of the right complement
          
          if (fn_r == 0) { sum_rght = 1;
          } else {
            
            // calculate initial recursion arguments
            r_level_init = fn_r;
            crt = (fn_l + 1) * c1 / c0 * n0 / n1 - tn_l;
            r_mxpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt + 1) : ceil(crt);
            r_mxpsr = std::max(0, r_mxpsr);
            r_start_init = r_minpsr + (fn_r - 1) + 1;
            r_stop_init = (fn_r + tn_r) - r_mxpsr;
            
            if (r_start_init <= r_stop_init & r_stop_init <= fn_r + tn_r) {
              
              // start recursion
              negComp_rec(memo, &sum_rght, r_level_init, r_start_init, r_stop_init, tn_r, tn_l, fn_r, fn_l, n0, n1, c0, c1, pi0);

            } else { sum_rght = 0; }
            
          }
          
          // overall number of permutations
          if (l == mid) {
            sum_tot = sum_tot + sum_lft * sum_rght;
          } else {
            sum_tot = sum_tot + 2 * sum_lft * sum_rght;
          }
      
        }
      }
    }
    
    perm[i] = sum_tot;
    
  }
  
  memo.clear(); // delete map
  
  return(wrap(perm));

}

