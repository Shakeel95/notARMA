##############################################
# functions templates for NOT-ARMA-1 and NOT-ARMA-2
# (for get.not.arma.bic.mod, it is in another file)
##############################################

get.predict <- function (not.arma.mod, test.i, ...) {
  UseMethod("get.predict", not.arma.mod)
}

get.mse <- function(not.arma.mod, test.d, test.i, ...) {
  pred <- get.predict(not.arma.mod, test.i = test.i, ...)
  return (mse.from.pred(pred = pred, test.d = test.d))
}

get.bic <- function(mod, not.mod, ...) {
  UseMethod("get.bic", mod)
}

# get the solution solution p
get.not.arma.sol.path <- function(not.mod, is.detrend) {
  mods <- list()
  all.cps <- get.all.not.ths.n.cps(not.mod)$all.cps
  for(i in 1:length(all.cps)) {
    cps <- sort(all.cps[[i]])
    if(is.detrend) {
      mods[[i]] <- get.not.arma.2.given.cps(not.mod = not.mod, cps = cps)
    } else {
      mods[[i]] <- get.not.arma.1.given.cps(not.mod = not.mod, cps = cps) 
    }
  }
  return (mods)
}