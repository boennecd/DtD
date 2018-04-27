#' @importFrom checkmate assert_numeric assert_number assert_choice
#' @export
BS_fit <- function(S, D, T., r, time, vol_start, method = c("kmv", "mle"),
                   tol = 1e-12, eps = 1e-8){
  #####
  # checks
  method <- method[1]
  assert_choice(method, c("kmv", "mle"))

  assert_numeric(S   , lower = 1e-16, finite = TRUE)
  assert_numeric(D   , lower = 1e-16, finite = TRUE)
  assert_numeric(T.  , lower = 1e-16, finite = TRUE)
  assert_numeric(r                  , finite = TRUE)
  assert_numeric(time               , finite = TRUE)
  stopifnot(!is.unsorted(time))

  if(missing(vol_start)){
    # use heuristic from Bharath et al. (2008)
    vol_start <- sd(diff(log(S))) * tail(S, 1) / (tail(S, 1) + tail(D, 1))

  } else
    assert_number(vol_start, lower = 1e-16, finite = TRUE)
  assert_number(tol      , lower = 1e-16, finite = TRUE)
  assert_number(eps      , lower = 1e-16, finite = TRUE)

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
