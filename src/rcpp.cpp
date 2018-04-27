#include "bs.h"

//' @export
// [[Rcpp::export]]
double BS_call(
    const double V, const double D, const double T, const double r,
    const double std){
  return BS_call_cpp(V, D, T, r, std);
}

// [[Rcpp::export]]
Rcpp::List BS_fit_cpp(
    const arma::vec &S, const arma::vec &D, const arma::vec &T,
    const arma::vec &r, const arma::vec &time, double vol_start,
    const std::string method,
    const double tol, const double eps){
  est_result res;

  if(method == "kmv"){
    res = kmv(S, D, T, r, time, vol_start, tol, eps);
    if(!res.success)
      Rcpp::stop("KMV method failed");

  } else
    Rcpp::stop("Method not implemented");

  return Rcpp::List::create(
    Rcpp::_["ests"] = Rcpp::NumericVector::create(
      Rcpp::_["mu"] = res.mu, Rcpp::_["vol"] = res.vol),
    Rcpp::_["n_iter"]  = res.n_iter,
    Rcpp::_["success"] = res.success);
}
