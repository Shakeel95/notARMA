get.not.arma.cv.mod <- function(not.mod, num.validate, is.detrend, cv.method, 
                                validate.fit = NULL,
                                return.validate.fit = FALSE) {
  target.i <- length(not.mod$x)
  arma.fit <- get.not.arma.bic.mod(not.mod, 
                                   is.detrend = is.detrend, 
                                   return.all.mods = TRUE)
  
  validate.fit <- update.validate.fit(d = not.mod$x, 
                                      target.not.mod = not.mod,
                                      num.validate = num.validate, 
                                      is.detrend = is.detrend, 
                                      contrast = not.mod$contrast, 
                                      validate.fit = validate.fit)
  
  th <- get.target.th(validate.fit = validate.fit, target.i, 
                      target.ths = arma.fit$all.ths,
                      num.validate = num.validate, cv.method = cv.method)
  mod <- list(mod = arma.fit$all.mods[[which(arma.fit$all.ths == th)]],
              th = th)
  if(return.validate.fit) {
    mod$validate.fit = validate.fit
  }
  return (mod)
}

update.validate.fit <- function(d, target.not.mod, num.validate, is.detrend, contrast, 
                                validate.fit = NULL) {
  validate.idx <- length(d) - c(num.validate:0) 
  
  all.idx <- union(validate.idx, validate.fit$idx)
  
  new.validate.fit <- list(is.detrend = is.detrend,
                           idx = all.idx)
  
  new.validate.fit$not.mods <- lapply(all.idx, function(validate.i) {
    i <- which(validate.i == validate.fit$idx)
    if(length(i)) validate.fit$not.mods[[i]] else NA 
  })
  
  new.validate.fit$arma.not.mse <- lapply(all.idx, function(validate.i) {
    i <- which(validate.i == validate.fit$idx)
    if(length(i)) validate.fit$arma.not.mse[[i]] else NA 
  })
  new.validate.fit$not.mod[[length(all.idx)]] <- target.not.mod
  new.validate.fit$arma.not.mse[[length(all.idx)]] <- NA
  
  for(i in 1:(length(all.idx)-1)) {
    validate.i <- all.idx[i]
    current.train.d <- d[1:(validate.i-1)]
    current.validate.d <- d[validate.i]
    
    if(is.na(new.validate.fit$not.mods[[i]])) {
      new.validate.fit$not.mods[[i]] <- get.not.mod(d = current.train.d, contrast = contrast)    
    }
    if(is.na(new.validate.fit$arma.not.mse[[i]])) {
      all.mods <- get.not.arma.bic.mod(not.mod = new.validate.fit$not.mods[[i]], 
                                       is.detrend = is.detrend, 
                                       return.all.mods = TRUE)$all.mods
      
      new.validate.fit$arma.not.mse[[i]] <- get.not.arma.sol.path.mse(not.mod = new.validate.fit$not.mods[[i]], 
                                                                      train.d = current.train.d, 
                                                                      test.d = current.validate.d, 
                                                                      not.arma.mods = all.mods)
    }
  }
  return (new.validate.fit)
}