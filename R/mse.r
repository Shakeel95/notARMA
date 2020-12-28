get.mse.n.cps.from.th.or.pos <- function(all.fits, target.th = NULL, target.pos.i = NULL, test.i) {
  if(is.null(target.pos.i) == is.null(target.th) || any(is.na(target.th))) {
    stop ("get.mse.from.th: must supply target.i or target.th")
  }
  
  extracted.fit <- extract.fits(all.fits, test.i)
  test.mse <- extracted.fit$mse[[1]]
  not.mod <- extracted.fit$not.mods[[1]]
  all.cps <- get.all.not.ths.n.cps(not.mod)$all.cps
  
  if(!is.null(target.th)) {
    ths <- get.all.not.ths.n.cps(not.mod)$all.ths
    target.pos.i <- which(ths==target.th)
  }
  
  return (list(mse = test.mse[target.pos.i],
               cps = all.cps[[target.pos.i]]))
}

mse.from.pred <- function(pred, test.d) {
  return (mean((pred-test.d)^2))
}

get.not.arma.sol.path.mse <- function(not.mod, train.d, test.d, not.arma.mods) {
  num.mod <- length(not.arma.mods)
  # calculate mse
  not.arma.mse <- rep(NA, num.mod)
  test.i <- get.test.idx.from.train.n.test.d(train.d, test.d)
  for(i in 1:num.mod) {
    not.arma.mse[i] <- get.mse(not.arma.mods[[i]], test.d = test.d, test.i = test.i, 
                               not.mod = not.mod)
  }
  return (not.arma.mse)
}

repartition.mse <- function(validate.fit, target.ths) {
  repartition.one.mse <- function(mse, ths, target.ths) {
    valid.i <- which(!is.na(mse))
    repartitioned.mse <- stats::approx(x = ths[valid.i], y = mse[valid.i], xout = target.ths)$y
    
    i <- which(is.na(repartitioned.mse))
    valid.i <- range(which(!is.na(repartitioned.mse)))
    
    large.ths.i <- i[i<valid.i[1]]
    repartitioned.mse[large.ths.i] <- mse[1]
    
    small.ths.i <- i[i>=valid.i[2]]
    repartitioned.mse[small.ths.i] <- tail(mse, 1)
    
    return (repartitioned.mse)
  }
  
  # each column is a thr
  return (do.call(rbind, lapply(1:length(validate.fit$idx), function(i) {
    not.mod <- validate.fit$not.mods[[i]]
    return (repartition.one.mse(mse = validate.fit$arma.not.mse[[i]], 
                                ths = get.all.not.ths.n.cps(not.mod)$all.ths, 
                                target.ths = target.ths))
  })))
}