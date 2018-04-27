#' @importFrom checkmate assert_numeric assert_number assert_choice
#' @export
BS_sim <- function(vol, mu, dt, V_0, D, r, T.){
  #####
  # checks
  assert_number(vol, lower = 1e-16, finite = TRUE)
  assert_number(mu                , finite = TRUE)
  assert_number(dt , lower = 1e-16, finite = TRUE)
  assert_number(V_0, lower = 1e-16, finite = TRUE)

  assert_numeric(D , lower = 1e-16, finite = TRUE)
  assert_numeric(r                , finite = TRUE)
  assert_numeric(T., lower = 1e-16, finite = TRUE)

  #####
  # get arguments for latter calls
  lens = c(length(D), length(T.), length(r), 1L)
  args <- .get_eq_length_args(
    lens = lens, D = D, T = T., r = r, vol = vol)

  #####
  # simulate and return
  time <- (1:max(lens) - 1L) * dt
  V <- V_0 * exp(
    (mu - vol^2/2) * time + cumsum(c(
      0, rnorm(length(time) - 1L, sd = vol * sqrt(dt)))))
  S <- with(args, mapply(BS_call, V, D, T = T, r, vol))

  with(args, data.frame(V, S, time, D, r, vol, mu, T))
}
