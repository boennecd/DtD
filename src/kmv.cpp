#include "bs.h"

inline bool almost_eq(double a, double b, double eps){
  return std::abs(a - b) / (std::abs(a) + 1e-8) < eps;
}

est_result kmv(
    const arma::vec &S, const arma::vec &D, const arma::vec &T,
    const arma::vec &r, const arma::vec &time, double vol_start,
    const double tol, const double eps){
  est_result out;

  // assume that all have equal length
  const arma::uword n = S.n_elem;

  // compute dts
  arma::vec dts(n - 1L);
  const double *time_i = time.begin();
  for(auto it = dts.begin(); it != dts.end(); ++it, ++time_i)
    *it = *(time_i + 1) - *time_i;

  // find solution
  double &mu = out.mu, &vol = out.vol;
  vol = vol_start;
  const signed int it_max = 1000L;
  arma::vec vol_vec(n);

  signed int &i = out.n_iter;
  for(i = 0; i < it_max; ++i){
    double vol_old = vol, mu_old = mu;
    vol_vec.fill(vol);
    arma::vec V = get_underlying_cpp(S, D, T, r, vol_vec, tol);

    // compute vol and mean
    double *dt = dts.begin();
    double log_prev = std::log(V[0]), log_new, s1 = 0, s2 = 0, ss = 0,
      log_return;
    for(auto V_i = V.begin() + 1; V_i != V.end(); ++V_i, ++dt){
      log_new = std::log(*V_i);
      log_return = log_new - log_prev;

      s1 += log_return              /           *dt;
      s2 += log_return              / std::sqrt(*dt);
      ss += log_return * log_return /           *dt;
      log_prev = log_new;
    }

    mu  =  s1 / (n - 1.);
    s2 /= (n - 1.);
    vol = std::sqrt(ss / (n - 1.) - s2 * s2);

    // check if converged
    if(i > 0L and almost_eq(mu_old, mu, eps) and
       almost_eq(vol_old, vol, eps)){
      out.success = true;
      ++i;
      break;
    }
  }

  return out;
}
