
test_that("calibration works", {
    disp <- weitrix_dispersions(simwei, ~1)
    cal <- weitrix_calibrate(simwei, disp)

    # Calibration should have adjusted weights so all dispersions equal 1
    new_disp <- weitrix_dispersions(cal, ~1)
    expect_equal(new_disp, rep(1, nrow(cal)))
})

