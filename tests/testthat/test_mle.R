context("Testing MLE method")

test_that("mle method gives previous output without supplying a starting value for the vol", {
  set.seed(83486778)
  sims <- BS_sim(
    vol = .1, mu = .05, dt = .1, V_0 = 100, T. = 1, D = rep(80, 20), r = .01)

  # with(sims, dput(BS_fit(S = S, D = D, T. = T, r = r, time = time, method = "mle")))
  with(sims,
       expect_equal(
         BS_fit(S = S, D = D, T. = T, r = r, time = time, method = "mle"),
         structure(list(
           ests = structure(c(0.0805269674101987, 0.101307053290089),
                            .Names = c("mu", "vol")),
           n_iter =  26L, success = TRUE),
           .Names = c("ests", "n_iter", "success"))))
})

test_that("mle method gives previous output when supplying a starting value for the vol", {
  set.seed(79156879)
  sims <- BS_sim(
    vol = .1, mu = .05, dt = .2, V_0 = 100, T. = 1, D = rep(80, 20), r = .01)

  # with(sims, dput(BS_fit(S = S, D = D, T. = T, r = r, time = time,
  #                        method = "mle", vol_start = 1)))
  with(sims,
       expect_equal(
         BS_fit(S = S, D = D, T. = T, r = r, time = time, method = "mle",
                vol_start = 1),
         structure(list(
           ests = structure(c(0.102790534402294, 0.118642486263979),
                            .Names = c("mu", "vol")),
           n_iter = 32L, success = TRUE),
           .Names = c("ests", "n_iter", "success"))))
})
