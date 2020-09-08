

test_that("least squares solver works", {
    set.seed(1234)
    n <- 10
    m <- 5
    X <- matrix(rnorm(n*m),nrow=n)
    y <- rnorm(n)
    w <- rnorm(n)^2
    present <- rep(TRUE, n)
    offset <- rep(0, n)

    func <- least_squares_func(X)

    # Should match weighted lm
    expect_equal(
        func(present, w, y, offset), as.vector(coef(lm(y ~ 0+X, weights=w))))

    # With new weights
    w <- rnorm(n)^2
    expect_equal(
        func(present, w, y, offset), as.vector(coef(lm(y ~ 0+X, weights=w))))

    # With reused weights
    expect_equal(
        func(present, w, y, offset), as.vector(coef(lm(y ~ 0+X, weights=w))))

    # With subset
    present[c(2,3)] <- FALSE
    expect_equal(func(present, w[present], y[present], offset), 
        as.vector(coef(lm(y[present] ~ 0+X[present,], weights=w[present]))))

    # With reused subset
    expect_equal(func(present, w[present], y[present], offset), 
        as.vector(coef(lm(y[present] ~ 0+X[present,], weights=w[present]))))
})
