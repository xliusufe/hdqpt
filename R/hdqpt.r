#' High dimensional quantile partial test
#'
#' @useDynLib hdqpt, .registration = TRUE
#' @export
#' @param x The design matrix, which is a n by p matrix.
#' @param y The outcome.
#' @param u The control factors, which is a n by q matrix.
#' @param tau The quantile.
#' @param coef The coefficient of control factors, the first element corresponding to intercept. If not provided, we will first estimate according to the user's choice of estimation method, and then test process.
#' @param decorrelate An indicator of whether do decorrelate test.
#' @param W The decorrelate matrix, should be a q by p matrix when provided. Necessary when decorrelate is TRUE. If not provided, we will estimate it using cv.glmnet function from glmnet.
#' @param method Method used to estimator the coefficient of control factors when it is not provided. "hqreg" means using "cv.hqreg" funcction from "hqreg" package and "rqPen" means using "rq.pen" function with "aLASSO" penalty.
#' @param seed seed for hqreg estimation process.
#' @importFrom glmnet cv.glmnet
#' @importFrom Matrix t
#' @importFrom stats pnorm quantile
#' @importFrom rqPen rq.pen qic.select
#' @importFrom hqreg cv.hqreg
#' @return A list.
#' \itemize{
#'     \item ts - Test statistic.
#'     \item pval - p value.
#'     \item tr_sigma - trace sigma square.
#'     \item coef - the coefficient of control factors.
#' }
#' @examples
#' set.seed(0)
#' n <- 300
#' p <- 710
#' x <- matrix(rnorm(n * p), n, p)
#' u <- matrix(rnorm(n * p), n, p)
#' alpha <- c(rep(1, 5), rep(0, p - 5))
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' y <-  u %*% alpha + x %*% beta + rnorm(n)
#' test_result <- hdqpt(x, u, y, tau = 0.5, decorrelate = FALSE, method = "hqreg")
#' print(test_result)
hdqpt <- function(x, u, y, tau = 0.5, coef = NULL, decorrelate = TRUE, W = NULL, method = c("hqreg", "rqPen"), seed = 1) {
    if (tau >= 1 || tau <= 0) stop('tau is out of range!')

    n  <- length(y)
    px <- ifelse(is.null(ncol(x)), 1, ncol(x))

    method <- match.arg(method)
    if (is.null(coef)) {
        if (method == "rqPen") {
            fit <- rqPen::qic.select(rqPen::rq.pen(u, y, tau = tau, penalty = "aLASSO"))
            coef <- as.vector(fit$coefficients)
            resids <- as.numeric(y - cbind(1, u) %*% coef)
        } else if (method == "hqreg") {
            fit <- hqreg::cv.hqreg(u, y, FUN = "hqreg_raw", seed = seed, message = FALSE)
            coef <- coef(fit)
            resids <- as.numeric(y - cbind(1, u) %*% coef)
        }
    } else {
        resids <- as.numeric(y - cbind(1, u) %*% coef)
    }
    resids <- resids - as.numeric(quantile(resids, probs = tau))
    resids <- ifelse(resids > 0, -tau, 1 - tau)

    if (decorrelate) {
        if (is.null(W)) {
            xhat <- NULL
            for (k in 1:px) {
                fitglm <- glmnet::cv.glmnet(u, x[, k])
                resids  <- x[, k] - u %*% coef(fitglm)[-1]
                xhat <- cbind(xhat, resids)
            }
        } else {
            xhat <- x - u %*% W
        }
    } else {
        xhat <- x
    }

    dims <- c(n, px)
    fit <- .Call(
        "_QPT_Test",
        as.numeric(t(xhat)),
        as.numeric(resids),
        as.integer(dims),
        as.numeric(tau)
    )

    result <- list(
        ts = fit$test,
        pvals = 1 - pnorm(fit$test),
        tr_sigma = fit$tr_sigma,
        coef = coef)

    return(result)
}

