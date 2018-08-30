#include "bs.h"
#include "utils.h"

inline bool almost_eq(double a, double b, double eps){
  double abs_a = std::abs(a);
  return (abs_a < eps) ?
    std::abs(a - b) < eps :
    std::abs(a - b) / (abs_a + 1e-8) < eps;
}

est_result est_iterative(
    const arma::vec &S, const arma::vec &D, const arma::vec &T,
    const arma::vec &r, const arma::vec &time, double vol_start,
    const double tol, const double eps){
  est_result out;

  // assume that all have equal length
  const arma::uword n = S.n_elem;

  // find solution
  const arma::vec dts = diff(time);
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
    const double *dt = dts.begin();
    double log_prev, log_new = std::log(V[0]), xbar = 0, log_return, i = 1.,
      sse = 0;

    for(auto V_i = V.begin() + 1; V_i != V.end(); ++V_i, ++dt, i += 1.){
      log_prev = log_new;
      log_new = std::log(*V_i);
      log_return = log_new - log_prev;

      /* stable methods as in
       *    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
       *
       *  r ~ N(\mu dt, dt \sigma^2)
       *  <=> r / dt ~ N(\mu, \sigma^2 / dt)
       *  <=> r / sqrt{dt} ~ N(\mu \sqrt{dt}, \sigma^2)
       *
       *  So first update the mean and then the sse
       */
      double xbar_old = xbar;
      xbar += (log_return  / *dt - xbar) / i;
      double t1 = std::sqrt(*dt), t2 = log_return / t1;
      sse += (t2 - xbar_old * t1) * (t2 - xbar * t1);
    }

    vol = std::sqrt(sse / (n - 1.)); /* recall this is the regular division by
                                        n as we only have n - 1 returns */
    mu = xbar + vol * vol / 2.;

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
