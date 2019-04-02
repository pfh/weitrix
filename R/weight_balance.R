

weight_balance <- function(weitrix, bio_weights=1, design=~1, p=0, verbose=TRUE, BPPARAM=getAutoBPPARAM()) {
    weitrix <- as_weitrix(weitrix)

    if (length(bio_weights) == 1)
        bio_weights <- rep(bio_weights, nrow(weitrix))
    assert_that(length(bio_weights) == nrow(weitrix))

    x <- weitrix_x(weitrix)
    tv_weights <- weitrix_weights(weitrix)
    present <- tv_weights>0

    # Find components assuming only biological noise
    bio_weitrix <- weitrix
    weitrix_weights(bio_weitrix) <- present * bio_weights
    comp <- weitrix_components(bio_weitrix, p=p, design=design, verbose=verbose, BPPARAM=BPPARAM)
    df_total <- comp$df_total
    df_residual <- comp$df_residual

    calc_var <- function(weights) {
        ss <- calc_weighted_ss(x, weights, comp$row, comp$col, BPPARAM=BPPARAM)
        ss / df_residual
    }

    score_weights <- function(param) {
        weights <- tv_weights/(1-param + param/bio_weights*tv_weights)
        result <- df_total*log(calc_var(weights)) - sum(log(weights[present]))
        cat(param,"gave",result,"\n")
        result
    }

    param <- optimize(score_weights, c(0, 1))$minimum

    balanced_weights <- tv_weights / (1-param+param/bio_weights*tv_weights)
    
    # Allow weights to be used directly as precisions
    residual_var <- calc_var(balanced_weights)
    metadata(weitrix)$weitrix$tech_scale <- residual_var*(1-param)
    metadata(weitrix)$weitrix$bio_scale <- residual_var*param
    weitrix_weights(weitrix) <- balanced_weights / residual_var

    weitrix
}
