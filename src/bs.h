#include "RcppArmadillo.h"
#include <array>

double BS_call_cpp(
    const double, const double, const double, const double, const double);

double BS_call_cpp_inv(
    const double, const double, const double, const double, const double,
    const double,
    double, double, double);

arma::vec get_underlying_cpp(
    const arma::vec&, const arma::vec&, const arma::vec&, const arma::vec&,
    const arma::vec&, const double);

struct est_result {
  double mu, vol;
  bool success;
  signed int n_iter;

  est_result():
    mu (std::numeric_limits<double>::quiet_NaN()),
    vol(std::numeric_limits<double>::quiet_NaN()),
    success(false), n_iter(0L) {}
  };

est_result kmv(
    const arma::vec&, const arma::vec&, const arma::vec&, const arma::vec &,
    const arma::vec&, double vol, const double, const double);
