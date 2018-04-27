#include "RcppArmadillo.h"
#include "bs.h"

//' @export
// [[Rcpp::export]]
double BS_call(
    const double V, const double D, const double T, const double r,
    const double std){
  return BS_call_cpp(V, D, T, r, std);
}
