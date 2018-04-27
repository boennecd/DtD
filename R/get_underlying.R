#' @export
get_underlying <- function(S, D, T., r, std, tol = 1e-12){
  lens <- c(length(S), length(D), length(T.), length(r), length(std))
  if(all(lens == 1))
    return(drop(get_underlying_cpp(
      S = S, D = D, T = T., r = r, std = std, tol = tol)))

  max_len <- max(lens)
  args <- list(S = S, D = D, T = T., r = r, std = std)
  if(max_len == 0 || any(!lens %in% c(1L, max_len)))
    stop("All arguments but ", sQuote("tol"), " should have the max length or",
         " length 1")

  not_max <- which(lens != max_len)
  args[not_max] <- lapply(args[not_max], rep, times = max_len)

  drop(get_underlying_cpp(
    S = args$S, D = args$D, T = args$T,
    r = args$r, std = args$std, tol = tol))
}
