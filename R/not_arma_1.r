################################
# NOT-ARMA-1 
################################

#' NOT-ARMA-1 
#' @export not_arma_1
#' @import not
#' @importFrom stats fitted lm predict
#' @param d data
#' @param contrast A type of the contrast function used in the NOT algorithm. Choice of "pcwsConstMean", "pcwsConstMeanHT", "pcwsLinContMean", "pcwsLinMean", "pcwsQuadMean", "pcwsConstMeanVar". See not::not for more details.
#' @param thr.method Threshold methods, "bic" for bic, "cv.1" for minimum cv mse, "cv.2" for 1se rule
#' @param num.validate optional, number of validation observations
#' @return a list with fitted NOT model, change points, BIC value and the threshold in NOT
#' @examples
#' \donttest{
#' d = get.s1.sim(n = 200, cp = 100, intercept = 0, jump.size = 1, sigma = 1, arma.mod = c(0.8, 0.2))$d
#' not_arma_1(d, contrast = "pcwsConstMean")
#' }
not_arma_1 <- function(d, contrast, thr.method = "bic", num.validate = round(length(d)*0.15)) {
  not.mod <- get.not.mod(d, contrast)
  if(thr.method == "bic") {
    return (get.not.arma.1.bic.mod(not.mod, return.all.mods = FALSE))
  } else {
      cv.method <- ifelse(thr.method == "cv.1",1,2)
      mod <- get.not.arma.1.cv.mod(not.mod = not.mod, cv.method = cv.method, 
                                   num.validate = num.validate)
    mod$not.mod <- not.mod
    return (mod)
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
#' d <- get.s1.sim(n = 200, cp = 100, intercept = 0, jump.size = 1, 
#' sigma = 1, arma.mod = c(0.8, 0.2))$d
#' mod <- not_arma_1(d[1:199], contrast = "pcwsConstMean")
#' pred <- predict(mod, 200)
#' }
predict.not.arma.1 <- function(object, test.i, ...) {
  return (get.predict(object$mod, test.i, ...))
}

# ======== helpers ============
get.not.arma.1.given.cps <- function(not.mod = NULL, cps, d = NULL, contrast = NULL) {
  if(!is.null(not.mod)) {
    d = not.mod$x
    contrast = not.mod$contrast
  } 
  if(is.null(d) || is.null(contrast)) {
    stop("get.not.arma.1.given.cps: provide not.mod or d and contrast")
  }
  
  end.pts <- get.left.n.right.end.pts(cps, length(d)) 
  mod <- lapply(1:length(end.pts$left), 
                function(i) {
                  one.mod <- get.approx.mod(x = d, start.i = end.pts$left[i], end.i = end.pts$right[i], 
                                            contrast = contrast)
                  if(is.null(get.fitted.approx.mod(one.mod))) {
                    stop("get.not.arma.1.given.cps: fail to get correct mod.")
                  }
                  return (list(mod = one.mod, end = end.pts$right[i]))
                })
  class(mod) <- "not.arma.1"
  return (mod)
}

get.predict.not.arma.1 <- function(not.arma.1.mod, test.i, ...) {
  segment.not.arma.1.mod <- not.arma.1.mod[[length(not.arma.1.mod)]]
  last.cp <- ifelse(length(not.arma.1.mod)>1, not.arma.1.mod[[length(not.arma.1.mod)-1]]$end, 0)
  return (get.predict.approx.mod(mod = segment.not.arma.1.mod$mod, test.i = test.i))
}

get.bic.not.arma.1 <- function(not.arma.1.mod, not.mod) {
  res <- unlist(lapply(not.arma.1.mod, 
                       function(segment.mod) get.residuals.approx.mod(segment.mod$mod)))
  p <- sum(sapply(not.arma.1.mod, function(segment.mod) sum(get.coef.approx.mod(segment.mod$mod)!= 0)))
  bic <- get.ic.from.mse(mse = mean(res^2), n = not.mod$n, p = p, 
                         m = length(not.arma.1.mod)-1, cp.penalty = 1, ic = "bic")  
  if(is.na(bic)) {
    stop("get.bic.not.arma.1: Inf bic")
  }
  return (bic)
}

get.not.arma.1.bic.mod <- function(not.mod, return.all.mods = FALSE) {
  mod <- get.not.arma.bic.mod(not.mod = not.mod, is.detrend = FALSE, return.all.mods = return.all.mods)
  class(mod) <- c("not.arma.1")
  return (mod)
}

get.not.arma.1.cv.mod <- function(not.mod, cv.method, num.validate = round(length(not.mod$x)*0.15)) {
  mod <- get.not.arma.cv.mod(not.mod = not.mod, 
                             num.validate = num.validate, is.detrend = FALSE, 
                             cv.method = cv.method, 
                             validate.fit = NULL,
                             return.validate.fit = FALSE)
  class(mod) <- c("not.arma.1")
  return (mod)
}

get.not.arma.1.sol.path <- function(not.mod) {
  return (get.not.arma.sol.path(not.mod = not.mod, is.detrend = FALSE))
}