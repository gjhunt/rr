center <- function(x) scale(x, scale = FALSE)

winsor <- function(x, fraction = 0.01) {
    if (length(fraction) != 1 || fraction < 0 || fraction > 0.5) {
        stop("bad value for 'fraction'")
    }
    lim <- stats::quantile(x, probs = c(fraction, 1 - fraction), na.rm = TRUE)
    x[x < lim[1]] <- lim[1]
    x[x > lim[2]] <- lim[2]
    x
}

trim <- function(x, alpha) {
    q <- stats::quantile(x, c(alpha, 1 - alpha))
    return(x[x > q[1] & x < q[2]])
}

gm_mean <- function(x) {
    exp(mean(log(x[x > 0])))
}
