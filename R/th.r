
get.target.th <- function(validate.fit, target.i, target.ths, num.validate, cv.method = 1) {
  validate.fit <- extract.fits(validate.fit, target.i-c(num.validate:1)) 
  repartitioned.mse <- repartition.mse(validate.fit = validate.fit, target.ths = target.ths)
  th <- get.cv.not.arma.given.validate.mse(target.ths = target.ths,
                                           validate.solution.path.mse = repartitioned.mse, 
                                           method = cv.method)
  return (th)
}

get.cv.not.arma.given.validate.mse <- function(target.ths, validate.solution.path.mse,
                                               method) {
  validate.mse.mean <- colMeans(validate.solution.path.mse)
  min.i <- min(which(validate.mse.mean == min(validate.mse.mean, na.rm = TRUE)))
  if(method == 1) {
    return (target.ths[min.i])
  }
  if(method == 2){
    validate.mse.se <- sapply(1:ncol(validate.solution.path.mse), 
                              function(i) stats::sd(validate.solution.path.mse[,i]))/sqrt(nrow(validate.solution.path.mse))
    th <- max(target.ths[validate.mse.mean <= (validate.mse.mean+validate.mse.se)[min.i] & !is.na(validate.mse.mean)])
    return (th)
  }
}

extract.fits <- function(validate.fit, idx) {
  match.i <- which(validate.fit$idx%in%idx)
  if(length(match.i) == length(idx)) {
    return (list(arma.not.mse = lapply(match.i, function(i) validate.fit$arma.not.mse[[i]]),
                 not.mods = lapply(match.i, function(i) validate.fit$not.mods[[i]]),
                 idx = idx,
                 is.detrend = validate.fit$is.detrend))  
  } else {
    stop("extract.fits: wrong idx")
  }
}