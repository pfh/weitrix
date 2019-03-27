

weight_balance <- function(errors, tech_weights, bio_weights) {
    n <- nrow(errors)
    m <- ncol(errors)
    
    if (length(bio_weights) == 1)
        bio_weights <- rep(bio_weights, n)
    stopifnot(length(bio_weights) == n)
        
    if (is.null(design)) 
        design <- cbind(rep(1,ncol(E)))
        
    stopifnot(nrow(design) == m)
    p <- ncol(design)
    
    # Ensure NAs have weight 0, then remove them
    tv_weights[is.na(E)] <- 0
    present <- tv_weights > 0
    n_present <- sum(present)
    df <- n_present - n*p
    E[!present] <- 0

    # Calculate Ordinary Least Squares residuals
    residuals2 <- map(seq_len(nrow(E)), function(i) {
        presenti <- present[i,]
        result <- rep(0, m)
        result[presenti] <- lm.fit(design[presenti,,drop=F], E[i,presenti])$residuals^2
        result
    })
    residuals2 <- do.call(rbind, residuals2)

    calc_var <- function(weights) {
        sum(residuals2*weights) / df
    }

    # Choose optimium weight for technical variance component
    #
    # variance model is: overall_variance * ((1-p)*tv + p*bv)
    # hence weight is: 1/((1-p)/tw+p/bw)
    #
    # This is based on Maximum Likelihood for the residuals (assumed normally distributed).
    # The ML overall_variance given the weights can be directly found and substituted in, yielding this optimization: 
    score_weights <- function(param) {
        weights <- tv_weights/(1-param + param/bio_weights*tv_weights)
        n_present*log(calc_var(weights)) - sum(log(weights[present]))
    }
    param <- optimize(score_weights, c(0, 1))$minimum

    elist$weights <- tv_weights / (1-param+param/bio_weights*tv_weights)
    
    # Allow weights to be used directly as precisions
    residual_var <- calc_var(elist$weights)
    elist$techvar <- residual_var*(1-param)
    elist$biovar <- residual_var*param
    elist$weights <- elist$weights / residual_var
    
    # Ensure weight 0 encoded as NA
    elist$E[elist$weights == 0] <- NA
    
    elist
}
