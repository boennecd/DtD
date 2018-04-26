#include <Rcpp.h>
#include "bs.h"

//' @export
// [[Rcpp::export]]
double BS_call(
    const double V, const double D, const double T, const double r,
    const double std){
  return BS_call_cpp(V, D, T, r, std);
}

//' @export
// [[Rcpp::export]]
double get_underlying(
    const double S, const double D, const double T, const double r,
    const double std, const double tol = 1e-12){
  double V_min = S;
  double V_max = 100 * S + D;
  double V_mid =  10 * S + D;
  return BS_call_cpp_inv(S, D, T, r, std, tol, V_min, V_max, V_mid);
}
