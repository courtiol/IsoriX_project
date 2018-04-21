test_that("main_fn", {
  cat("\n Testing main fn... \n")
  GNIPDataDEagg <- prepsources(data = GNIPDataDE)
  GermanFit <- isofit(data = GNIPDataDEagg, mean_model_fix = list(elev = TRUE, lat_abs = TRUE), verbose = FALSE)
  suppressWarnings(GermanScape <- isoscape(raster = ElevRasterDE, isofit = GermanFit, verbose = FALSE))
  expect_equal(GermanScape$isoscapes$mean[50, 40][[1]], expected = -55.583527288361878504)
  expect_equal(GermanScape$isoscapes$mean_predVar[50, 40][[1]], expected = 140.10497493914493816)
  expect_equal(GermanScape$isoscapes$mean_residVar[50, 40][[1]], expected = 293.97314939810542)
})

