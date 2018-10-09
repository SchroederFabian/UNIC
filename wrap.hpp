#include <RcppCommon.h>
#include <gmpxx.h>
//#include <gmp.h>

namespace Rcpp {

    template <> SEXP wrap(const std::vector<mpz_class> &);
}
 
 
#include <Rcpp.h>
//#include <gmpxx.h>
//#include <gmp.h>
 
namespace Rcpp {

  // [[Rcpp::export]]
  template <> SEXP wrap(const std::vector<mpz_class> &x) {
    
    int n = x.size();
    std::vector <std::string> y(n);
    
    for (int i = 0; i < n; i++) {
      y[i] = x[i].get_str();
    }
    Function bigz("as.bigz", "gmp");
    
    return bigz(y);
    
   }
}

