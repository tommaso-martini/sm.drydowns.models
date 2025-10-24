test_that("muSEC mono runs and respects bounds", {
  set.seed(1)
  theta0 <- 0.30
  dz     <- 0.10
  ET     <- rep(3, 48)  # mm/h for 2 days
  
  res <- musec_simulate(
    theta_init = theta0,
    structure  = list(dz = dz),
    Ks=0.5, n=3, theta_s=0.40, theta_r=0.05, theta_fc=0.25, theta_wp=0.10, theta_star=0.30,
    ET=ET, timestep="hourly", exp_drain=TRUE
  )
  
  expect_equal(res$mode, "mono")
  th <- res$series$theta
  expect_true(all(th <= theta0 + 1e-12))
  expect_true(all(th >= 0.05 - 1e-12))
})

test_that("muSEC multi runs and produces matrices", {
  set.seed(1)
  structure <- compute_dz_depths(Z_t = 0.6, dz_top = 0.05, dz_rest = 0.10, n_top_layers = 2)
  depths <- structure$depths; dz <- structure$dz
  theta_init <- infer_profile_from_averages(depths, dz, c("0.15"=0.28, "0.60"=0.30))
  ET <- rep(3, 24)
  
  res <- musec_simulate(
    theta_init = theta_init,
    structure  = list(depths = depths, dz = dz),
    Ks=0.5, n=3, theta_s=0.40, theta_r=0.05, theta_fc=0.25, theta_wp=0.10, theta_star=0.30,
    ET=ET, timestep="hourly", exp_drain=TRUE
  )
  
  expect_equal(res$mode, "multi")
  expect_true(is.matrix(res$series$theta_mm))
  expect_equal(nrow(res$series$theta_mm), length(theta_init))
  expect_equal(ncol(res$series$theta_mm), length(ET))
})