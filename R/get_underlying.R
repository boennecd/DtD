#' @export
get_underlying <- function(S, D, T., r, vol, tol = 1e-12){
  lens <- c(length(S), length(D), length(T.), length(r), length(vol))
  if(all(lens == 1)) # return quickly
    return(drop(get_underlying_cpp(
      S = S, D = D, T = T., r = r, vol = vol, tol = tol)))

  # make sure all args have same length and call cpp code
  args <- .get_eq_length_args(
    lens = lens, S = S, D = D, T = T., r = r, vol = vol)
  with.default(args, drop(get_underlying_cpp(
    S = S, D = D, T = T, r = r, vol = vol, tol)))
}
