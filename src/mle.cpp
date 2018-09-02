#include "optim.h"
#include "bs.h"
#include "utils.h"
#include <math.h>

class log_like {
  const arma::uword n;
  const arma::vec &S, &D, &T, &r, dts, log_D, log_dts, sqrt_ts;
  arma::vec vol_vec;

  const double tol;

public:
  log_like(
    const arma::vec &S, const arma::vec &D, const arma::vec &T,
    const arma::vec &r, const arma::vec time, const double tol):
  n(S.n_elem), S(S), D(D), T(T), r(r), dts(diff(time)), log_D(arma::log(D)),
  log_dts(arma::log(dts)), sqrt_ts(arma::sqrt(T)), vol_vec(n), tol(tol) {}

  double compute(const double mu, const double vol){
    // find underlying
    vol_vec.fill(vol);
    const arma::vec V = get_underlying_cpp(S, D, T, r, vol_vec, tol);

    // compute log likelihood
    double log_like_t1 = 0, log_like_t2 = 0, log_prev,
      log_new = std::log(V[0]);
    const double *dt = dts.begin(), *log_dt = log_dts.begin(),
      vol_sq = vol * vol, *log_D_i = log_D.begin() + 1, *r_i = r.begin() + 1,
      *t_i = T.begin() + 1, *sqrt_t = sqrt_ts.begin() + 1;
    for(auto V_i = V.begin() + 1;
        V_i != V.end();
        ++V_i, ++dt, ++log_dt, ++log_D_i, ++r_i, ++t_i, ++sqrt_t){
      log_prev = log_new;
      log_new = std::log(*V_i);

      double log_mean = (mu - vol_sq / 2) * *dt;
      double err = log_new - log_prev - log_mean;

      // Notice the t_i instead of dt
      double log_deriv = R::pnorm(
        (log_new - *log_D_i + (*r_i + vol_sq / 2) * *t_i) /
          (*sqrt_t * vol),
        0, 1, 1, 1);

      // the log(dt) is from the normalization constant
      log_like_t1 -= (err * err/ (vol_sq * *dt) + *log_dt);
      log_like_t2 -= (log_new + log_deriv);
    }

    return - ((double)n - 1.) * std::log(2. * M_PI * vol_sq) / 2. +
      log_like_t1 / 2 + log_like_t2;
  }
};

/*
 * Use Nelder-Mead method from `optim`. See
 *    r-source/src/library/stats/R/optim.R
 *    r-source/src/library/stats/src/optim.c#L244
 *    r-source/src/appl/optim.c
 */
static double optimfunc(int n, double *p, void *ex)
{
  log_like* ll = (log_like*) ex;

  /* -1 as minimization is the default */
  return -ll->compute(p[0], std::exp(p[1]));
}

est_result mle(
    const arma::vec &S, const arma::vec &D, const arma::vec &T,
    const arma::vec &r, const arma::vec &time, double vol_start,
    const double tol, const double eps){
  // assign arguments
  double dpar[2] = { 0.0 /* mu */, log(vol_start) /* log to keep pos */ },
         *opar = new double[2], val;
  int fail, fncount;
  const double alpha = 1;
  const double beta  = 0.5;
  const double gamm  = 2.0;

  est_result out;
  log_like ll(S, D, T, r, time, tol);
  optim(2 /* npar */, dpar, opar, &val, optimfunc, &fail, -1e8 /* abstol */,
        eps /* reltol */, (void *) &ll, alpha, beta, gamm, 0L /* trace */,
        &fncount, 10000L /* maxit */);

  out.mu  = opar[0];
  out.vol = std::exp(opar[1]);
  delete[] opar;
  out.n_iter = fncount;
  out.success = fail == 0L;

  return out;
}

// [[Rcpp::export]]
double merton_ll_cpp(
    const arma::vec &S, const arma::vec &D, const arma::vec &T,
    const arma::vec &r, const arma::vec &time,
    const double vol, const double mu, const double tol){
  return log_like(S, D, T, r, time, tol).compute(mu, vol);
}
