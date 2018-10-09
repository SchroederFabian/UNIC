/* These functions calculate the number of favorable permutations of the null distribution
 * of the eic classifier for a positive complement for a given set of parameters fp, fn, 
 * n0, n1, c0, c1, pi0. posComp calculates the starting values and calls the recursive 
 * function posComp_rec. 
 */

#include <wrap.hpp>
#include <map>
#include <gmpxx.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

void posComp_rec(std::map <std::vector<int>, mpz_class>& memo, mpz_class* sums, int level, int start, int stop, int tp_h, int tp_o, int fp_h, int fp_o, int n0, int n1, double c0, double c1, double pi0) {
  
  int arr[7] = {level, start, stop, tp_h, tp_o, fp_h, fp_o}; // map key
  std::vector<int> itm(arr, arr+7);
  
  if (memo.find(itm) != memo.end()) { // map look up
    
    *sums = *sums + memo.at(itm);
    
  } else { // if key does not exist in map, calculate the corresponding value
    
    if (level == 1) { // end of recursion
      
      *sums = *sums + stop - start + 1;
      
    } else { // recursion step
      
      // calculate new recursion arguments
      int level_new = level - 1;
      double crt = level_new * c0 / c1 * n1 / n0;
      int mnpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
      crt = std::max((fp_o + (fp_h - level_new + 1)) * c0 / c1 * n1 / n0 - tp_o, 0.0); 
      int mxpsr = fabs(round(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
      int start_new = mnpsr + level_new;
      int stop_new = (fp_h + tp_h) - (mxpsr + (fp_h - level_new));
      
      
      mpz_class bfr(*sums);
      
      for (int i = start; i <= stop; i = i + 1) {
        posComp_rec(memo, sums, level_new, start_new, std::min(stop_new, i - 1), tp_h, tp_o, fp_h, fp_o, n0, n1, c0, c1, pi0);
      }
      
      mpz_class aftr(*sums);
      memo.insert( std::make_pair(itm, aftr - bfr)); // save pair of (key, value) into map
      
    }
  }
}


// [[Rcpp::export]]
SEXP posComp(Rcpp::NumericMatrix clst, int n0, int n1, double c0, double c1, double pi0) {

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
    
    crt = static_cast<double>(n1 - fn + fp - 1) / 2;
    int mx = static_cast<int>(ceil(crt));
    mid = ((n1 - fn + fp - 1) % 2 == 0) ? -1 : mx;
    
    for (int l = 1; l <= mx; l = l + 1) {
      for (int fp_l = 0; fp_l <= fp; fp_l = fp_l + 1) {
        
        int tp_l = l - fp_l; // true negatives on the left side
        int tp_r = (n1 - fn) - tp_l; // true negatives on the right side
        int fp_r = fp - fp_l; // false negatives on the right side

        sum_lft = 0;
        sum_rght = 0;
        
        // set l_minpsr and r_minpsr        
        if (fp_l == 0) { l_minpsr = -1;
        } else {
          crt = fp_l * c0 / c1 * n1 / n0;
          l_minpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
        }

        if (fp_r == 0) { r_minpsr = -1;
        } else {
          crt = fp_r * c0 / c1 * n1 / n0;
          r_minpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
        }
        
        if (tp_l >= l_minpsr & tp_r >= r_minpsr) { 
          
          // permutations of the left complement
          if (fp_l == 0) { sum_lft = 1;
          } else {
            
            
            // calculate initial recursion arguments
            l_level_init = fp_l;
            crt = (fp_r + 1) * c0 / c1 * n1 / n0 - tp_r; 
            l_mxpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt); 
            l_mxpsr = std::max(0, l_mxpsr);
            l_start_init = l_minpsr + (fp_l - 1) + 1;
            l_stop_init = (fp_l + tp_l) - l_mxpsr;
            
            if (l_start_init <= l_stop_init & l_stop_init <= fp_l + tp_l) {
              
              // start recursion
              posComp_rec(memo, &sum_lft, l_level_init, l_start_init, l_stop_init, tp_l, tp_r, fp_l, fp_r, n0, n1, c0, c1, pi0);

            } else { sum_lft = 0; }
          }
          
        // permutations of the right complement 
          
        if ( fp_r == 0) { sum_rght = 1;
        } else {
          
          // calculate initial recursion arguments
          r_level_init = fp_r;
          crt = (fp_l + 1) * c0 / c1 * n1 / n0 - tp_l;
          r_mxpsr = fabs(ceil(crt) - crt) < 1e-6 ? round(crt) : ceil(crt);
          r_mxpsr = std::max(r_mxpsr, 0);
          r_start_init = r_minpsr + (fp_r - 1) + 1;
          r_stop_init = (fp_r + tp_r) - r_mxpsr;
        
          if (r_start_init <= r_stop_init & r_stop_init <= fp_r + tp_r) {
            
            // start recursion
            posComp_rec(memo, &sum_rght, r_level_init, r_start_init, r_stop_init, tp_r, tp_l, fp_r, fp_l, n0, n1, c0, c1, pi0);

          } else { sum_rght = 0; }
          
        }
        
        // overall number of permutations equals their product
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
