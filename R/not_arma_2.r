################################
# NOT-ARMA-2
################################

#' NOT-ARMA-2
#' @export not.arma.2
#' @import not
#' @importFrom stats fitted lm predict
#' @param d data
#' @param contrast A type of the contrast function used in the NOT algorithm. Choice of "pcwsConstMean", "pcwsConstMeanHT", "pcwsLinContMean", "pcwsLinMean", "pcwsQuadMean", "pcwsConstMeanVar". See not::not for more details.
#' @param thr.method Threshold method. Currently only BIC is implemented
#' @return a list with fitted NOT model, change points, BIC value and the threshold in NOT
#' @examples
#' \donttest{
#' d = get.s1.sim(n = 200, cp = 100, intercept = 0, jump.size = 1, sigma = 1, arma.mod = c(0.8, 0.2))$d
#' not.arma.2(d, contrast = "pcwsConstMean")
#' }
not.arma.2 <- function(d, contrast, thr.method = "bic") {
  not.mod <- get.not.mod(d, contrast)
  if(thr.method == "bic") {
    mod <- get.not.arma.2.bic.mod(not.mod, return.all.mods = FALSE)
    mod$not.mod <- not.mod
    return (mod)
  } else {
    return (NULL)
  }
}

#' one step ahead prediction form fitted NOT-ARMA-1 model
#' @export
#' @param object fitted model from not.arma.1()
#' @param test.i time point to predict
#' @param ... additional arguments
#' @return one step ahead value
#' @examples
#' \donttest{
#' d <- get.s1.sim(n = 200, cp = 100, intercept = 0, jump.size = 1, sigma = 1, arma.mod = c(0.8, 0.2))$d
#' mod <- not.arma.2(d[1:199], contrast = "pcwsConstMean")
#' pred <- predict(mod, 200)
#' }
predict.not.arma.2.bic <- function(object, test.i, ...) {
  return (get.predict(object$mod, test.i, object$not, ...))
}

get.predict.not.arma.2 <- function(not.arma.2.mod, test.i, not.mod, ...) {

  not.pred <- get.not.signal.predict(not.mod, cps = not.arma.2.mod$cps, h = length(test.i))
  arma.pred <- get.predict.approx.mod(not.arma.2.mod$arma.mod, test.i = test.i)
  return (not.pred + arma.pred)
}

# ======== helpers ============

get.not.arma.2.given.cps <- function(not.mod, cps) {
  not.fitted <- predict(not.mod, cpt = cps) # from NOT
  detrend.x <- not.mod$x-not.fitted
  arma.mod <- get.approx.mod(x = detrend.x, contrast = not.mod$contrast)
  if(!is.null(get.fitted.approx.mod(arma.mod)) && length(get.fitted.approx.mod(arma.mod)) == length(not.fitted)) {
    arma.fit <- get.fitted.approx.mod(arma.mod)
  } else {
    arma.fit <- detrend.x - get.residuals.approx.mod(arma.mod)
  }
  result <- list(arma.mod = arma.mod,
                 fitted = not.fitted + arma.fit,
                 cps = cps)
  class(result) <- "not.arma.2"
  return (result)
}

get.bic.not.arma.2 <- function(not.arma.2.mod, not.mod) {
  res <- not.mod$x - not.arma.2.mod$fitted
  p <- sum(get.coef.approx.mod(not.arma.2.mod$arma.mod) != 0)
  cps <- not.arma.2.mod$cps
  num.cp <- ifelse(length(cps)== 1 && cps == 0, 0, length(cps))
  if(not.mod$contrast == "pcwsConstMean") {
    p <- p + (num.cp+1) # (num.cp+1) slopes
  } else if(not.mod$contrast == "pcwsLinContMean") {
    p <- p + 1 + (num.cp+1) # 1 intercepts and (num.cp+1) slopes
  } else if(not.mod$contrast == "pcwsLinMean") {
    p <- p + (num.cp+1)*2 # (num.cp+1) intercepts and slopes
  }
  bic <- get.ic.from.mse(mse = mean(res^2), n = not.mod$n, p = p,
                         m = num.cp, cp.penalty = 1, ic = "bic")
  if(is.na(bic)) {
    stop("get.bic.not.arma.2: Inf bic")
  }
  return (bic)
}

get.not.arma.2.bic.mod <- function(not.mod, return.all.mods = FALSE) {
  mod <- get.not.arma.bic.mod(not.mod = not.mod, is.detrend = TRUE, 
                              return.all.mods = return.all.mods)
  class(mod) <- c("not.arma.2.bic")
  return (mod)
}

get.not.arma.2.sol.path <- function(not.mod) {
  return (get.not.arma.sol.path(not.mod = not.mod, is.detrend = TRUE))
}

# =========== helper ============

# this predict trend (without arma part)
get.not.signal.predict <- function(not.mod, cps, h = 1) {
  fitted.val <- predict(not.mod, cpt = cps) # from NOT
  n <- length(fitted.val)
  cps <- cps[which(cps != n)]
  start.i <- ifelse(is.null(cps) || length(cps) == 0 || (length(cps) == 1 & cps == 0), 1, (max(cps)+1))
  current.trend <- fitted.val[start.i:n]
  slope <- mean(diff(current.trend))
  return (fitted.val[n]+cumsum(rep(slope, h)))
}