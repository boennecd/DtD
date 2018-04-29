#' @title  Fit Black-Scholes Parameters
#'
#' @description
#' Function to estimate the volatility, \eqn{\sigma}, and drift, \eqn{\mu}. See [TODO: insert vignette link]
#' for details. All vectors with length greater than one needs to have the same
#' length. The Nelder-Mead method from \code{\link{optim}} is used when
#' \code{method = "mle"}.
#'
#' @param S numeric vector with observed stock prices.
#' @param D numeric vector or scalar with debt due in \code{T.}.
#' @param T. numeric vector or scalar with time to maturity.
#' @param time numeric vector with the observation times.
#' @param vol_start numeric scalar with starting value for \eqn{\sigma}.
#' @param method string to specifiy which estimation method to use.
#' @param tol numeric scalar with tolerance in \code{\link{get_underlying}}.
#' @param eps convergence threshold.
#'
#' @return
#' A list with the following components
#' \item{ests}{estimates of \eqn{\sigma}, and drift, \eqn{\mu}.}
#' \item{n_iter}{number of iterations when \code{method = "iterative"}
#' and number of log likelihood evaluations when \code{method = "mle"}.}
#' \item{success}{logical for whether the estimation method converged.}
#'
#' @examples
#' library(DtD)
#' set.seed(83486778)
#' sims <- BS_sim(
#'   vol = .1, mu = .05, dt = .1, V_0 = 100, T. = 1, D = rep(80, 20), r = .01)
#'
#' with(sims,
#'      BS_fit(S = S, D = D, T. = T, r = r, time = time, method = "mle"))
#'
#' @importFrom checkmate assert_number assert_choice
#' @export
BS_fit <- function(S, D, T., r, time, vol_start, method = c("iterative", "mle"),
                   tol = 1e-12, eps = 1e-8){
  #####
  # checks
  method <- method[1]
  assert_choice(method, c("iterative", "mle"))
  .check_args(S = S, D = D, T. = T., r = r, time = time, tol = tol, eps = eps)

  if(missing(vol_start)){
    # use heuristic from Bharath et al. (2008)
    vol_start <- sd(diff(log(S))) * tail(S, 1) / (tail(S, 1) + tail(D, 1))

  } else
    assert_number(vol_start, lower = 1e-16, finite = TRUE)

  # some more checks when we have the lengths
  lens <- c(length(S), length(D), length(T.), length(r), length(time))
  max_len = max(max(lens))
  stopifnot(length(time) == max_len)
  stopifnot(max_len > 2L)

  #####
  # call cpp code and return
  args <- .get_eq_length_args(
    lens = lens, S = S, D = D, T = T., r = r, time = time)
  with.default(args, BS_fit_cpp(S, D, T, r, time, vol_start, method, tol, eps))
}
