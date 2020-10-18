##############################################
# ARIMA wrappers
# getter functions are implemented to handle possible
# missing data
##############################################
get.approx.mod <- function(x, p = 1, q = 1, start.i = 1, end.i = length(x), 
                           contrast, seed = NULL) {
  if(end.i < start.i | start.i < 1 | end.i > length(x)) {
    stop("get.approx.mod: wrong start.i or end.i")
  }
  if(!is.null(seed)) {
    set.seed(seed)
  }
  mod <- get.arima.mod.given.order(x = x, start.i = start.i, end.i = end.i, 
                                   contrast = contrast, seed = seed)

  if(is.null(stats::fitted(mod))) {
    stop ("get.approx.mod: fit error, fitted is null.")
  }
  wrap.up.mod <- list(mod = mod,
                      fitted = stats::fitted(mod),
                      xreg = mod$xreg)
  class(wrap.up.mod) <- "approx.mod"
  return (wrap.up.mod)
}

get.fitted.approx.mod <- function(mod) {
  if(is.null(mod$fitted)) {
    stop ("get.fitted.approx.mod: no fitted")
  }
  return (mod$fitted)
}

get.residuals.approx.mod <- function(mod) {
  if(is.null(mod$mod$residuals)) {
    stop ("get.coef.approx.mod: no residuals")
  }
  return (mod$mod$residuals)
}

get.coef.approx.mod <- function(mod) {
  if(is.null(mod$mod$coef)) {
    stop ("get.coef.approx.mod: no coef")
  }
  return (mod$mod$coef)
}

get.predict.approx.mod <- function(mod, test.i) {
  if(class(mod$mod) == "lm") {
    return (stats::predict(mod$mod, newdata = data.frame(xreg = test.i)))
  } else if(class(mod$mod) == "Arima") {
    xreg <- mod$xreg
    if(!is.null(xreg)) {
      return (as.numeric(stats::predict(mod$mod, newxreg = test.i)$pred))
    } else {
      return (as.numeric(stats::predict(mod$mod, n.ahead = length(test.i))$pred))
    }
  } else {
    stop(paste("get.predict.approx.mod: wrong model class", class(mod$mod)))
    return (NULL)
  }
}

get.arima.mod.given.order.helper <- function(x, order, xreg, contrast) {
  is.const.mean <- contrast == "pcwsConstMean" || contrast == "pcwsConstMeanVar"
  for(optim.method in c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN","Brent")) {
    if(is.const.mean) {
      one.mod <- tryCatch(stats::arima(x = x, order = order, optim.method = optim.method, 
                                       include.mean = TRUE),
                          error = function(e) NULL)
    } else{
      one.mod <- tryCatch(stats::arima(x = x, order = order, xreg = xreg, 
                                       optim.method = optim.method, include.mean = TRUE),
                          error = function(e) NULL)
    }
    if(!is.null(one.mod)) {
      if(!is.const.mean) { 
        one.mod$xreg <- xreg # make sure xreg is available
      }
      if(is.null(stats::fitted(one.mod))) {
        one.mod$fitted <- x - one.mod$residuals # make sure fitted is available
      }
      return (one.mod)
    }
  }
  return (NULL) # return NULL if cannot optimize...    
}

get.arima.mod.given.order <- function(x, p = 1, q = 1, d = 0, start.i = 1, end.i = length(x), 
                                      contrast, seed = NULL) {
  x.segment <- x[start.i:end.i]
  current.p <- p
  current.q <- q
  
  # try decreasing order of ARIMA till a valid model is returned
  while(current.p >= 0 & current.q >= 0) {
    mod <- get.arima.mod.given.order.helper(x = x.segment, order = c(current.p, d, current.q), 
                        xreg = start.i:end.i, contrast = contrast)
    if(!is.null(mod)) {
      return (mod)
    }
    
    if(current.p > 0) {
      new.order <- c((current.p-1), d, current.q)
      print("get.arima.mod.given.order: need to try other orders:")
      print(new.order)
      mod <- get.arima.mod.given.order.helper(x = x.segment, order = new.order, 
                          xreg = start.i:end.i, contrast = contrast)
      if(!is.null(mod)) return (mod)
    }
    if(current.q > 0) {
      new.order <- c(current.p, d, (current.q-1))
      print("get.arima.mod.given.order: need to try other orders:")
      print(new.order)
      mod <- get.arima.mod.given.order.helper(x = x.segment, order = new.order, 
                         xreg = start.i:end.i, contrast = contrast)
      if(!is.null(mod)) return (mod)
    }
    current.p <- current.p - 1
    current.q <- current.q - 1
  }
  
  # use p,q = 0 if all other possible order failed
  print("get.arima.mod.given.order: no arma model found, fit lm...")
  x.segment.d <- data.frame(y = x.segment, xreg = start.i:end.i)
  if(contrast == "pcwsConstMean" || contrast == "pcwsConstMeanVar") {
    mod <- lm(y~1, x.segment.d)
  } else {
    mod <- lm(y~xreg, x.segment.d)
  }
  return (mod)
}
