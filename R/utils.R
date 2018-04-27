#' @importFrom checkmate assert_integer
.get_eq_length_args <- function(lens, ...){
  args <- list(...)
  assert_integer(lens, lower = 1L, len = length(args))

  max_len <- max(lens)
  if(max_len == 0 || any(!lens %in% c(1L, max_len)))
    stop("Argument lengths are not valid")

  not_max <- which(lens != max_len)
  args[not_max] <- lapply(args[not_max], rep, times = max_len)
  args
}


