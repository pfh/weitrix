

# 

test_that("weitrix_confects behaves", {
    set.seed(1)
    for(i in seq_len(1)) {
        n <- 50
        m <- 10
        p <- 3
        design <- matrix(rnorm(m*p),nrow=m)   # random design matrix
        beta <- matrix(rnorm(n*p),nrow=n)     # random coefficients
        scaling <- (rchisq(n, 10)/10)^-0.5    # prior df 10
        y <- t(design %*% t(beta)) +          
            matrix(rnorm(n*m),nrow=n)*scaling # random residuals

        res1 <- weitrix_confects(y, design, c(0,0,1), step=0.001, full=TRUE)

        # ==== Results should exactly match limma ====
        library(limma)
        library(topconfects)
        fit <- lmFit(y, design)
        res2 <- limma_confects(fit, coef=3, step=0.001, full=TRUE)

        t1 <- res1$table[res1$table$index,]
        t2 <- res2$table[res2$table$index,]
        expect_equal( t1$effect, t2$effect )
        expect_equal( t1$confect, t2$confect )
        expect_equal( t1$fdr_zero, t2$fdr_zero )

        # ==== sd confects should be invariant to equivalent contrasts ====
        res3 <- weitrix_confects(y, design, 
            cbind(foo=c(0,1,0),bar=c(0,0,1)), full=TRUE)

        res4 <- weitrix_confects(y, design, 
            cbind(foo=c(0,1,1),bar=c(0,-1,1)), full=TRUE)

        t1 <- res3$table[res1$table$index,]
        t2 <- res4$table[res2$table$index,]
        expect_equal( t1$effect, t2$effect )
        expect_equal( t1$confect, t2$confect )
        expect_equal( t1$fdr_zero, t2$fdr_zero )

        # ==== Cohen's f confects should be invariant to equivalent contrasts ==
        res5 <- weitrix_confects(y, design, 
            cbind(foo=c(0,1,0),bar=c(0,0,1)), effect="cohen_f", full=TRUE)

        res6 <- weitrix_confects(y, design, 
            cbind(foo=c(0,1,1),bar=c(0,-1,1)), effect="cohen_f", full=TRUE)

        t1 <- res5$table[res1$table$index,]
        t2 <- res6$table[res2$table$index,]
        expect_equal( t1$effect, t2$effect )
        expect_equal( t1$confect, t2$confect )
        expect_equal( t1$fdr_zero, t2$fdr_zero )

        # TODO: further checking
    }
})



