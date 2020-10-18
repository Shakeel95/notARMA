#' simulate piecewise constant mean with ARMA noise
#' @export get.s1.sim
#' @importFrom stats arima.sim rnorm
#' @importFrom utils tail
#' @param n length of time series
#' @param cp change points
#' @param intercept starting mean
#' @param jump.size change(s) in mean
#' @param sigma (optional) if ARMA noise is not used
#' @param arma.mod (optional) ARMA structure of the noise
#' @return simulated data
#' @examples
#' \donttest{
#' d = get.s1.sim(n = 200, cp = 100, intercept = 0, jump.size = 1, sigma = 1, arma.mod = c(0.8, 0.2))
#' }
get.s1.sim <- function(n, cp, intercept, jump.size, sigma = NULL, arma.mod = NULL) { 
  # jump: yes, change in slope: no (make use of s3)
  return (get.s3.sim(n = n, cp = cp, 
                     start.slope = 0, change.slope = rep(0, length(cp)), 
                     intercept = intercept, jump.size = jump.size,
                     sigma = sigma, arma.mod = arma.mod))
}

# jump: no, change in slope: yes (make use of s3)
get.s2.sim <- function(n, cp, change.slope, start.slope, intercept, sigma = NULL, arma.mod = NULL) { 
  return (get.s3.sim(n = n, cp = cp, change.slope = change.slope,
                     start.slope = start.slope, intercept = intercept,
                     jump.size = rep(0, length(cp)),
                     sigma = sigma, arma.mod = arma.mod))
}

# jump: yes, change in slope: yes, change in vol: no (make use of helper)
get.s3.sim <- function(n, cp, start.slope, change.slope, intercept, jump.size, sigma = NULL, arma.mod = NULL) { 
  return (get.not.sim.helper(n = n, cp = cp, 
                             start.slope = start.slope, change.slope = change.slope,
                             intercept = intercept, jump.size = jump.size,
                             init.sigma = sigma, arma.mod = arma.mod))
}

# jump: yes, change in slope: no, change in vol: yes (make use of helper)
get.s4.sim <- function(n, cp, intercept, jump.size, sigma, jump.sigma) { 
  return (get.not.sim.helper(n = n, cp = cp, 
                             change.slope = rep(0, length(cp)), start.slope = 0, 
                             intercept = intercept, jump.size = jump.size,
                             init.sigma = sigma, jump.sigma = jump.sigma))
}

# ============ arma ==========
my.arma.sim <- function(mod, n, prev.x = NULL, prev.e = NULL, burn.period = 0) {
  ar <- mod$ar
  ma <- mod$ma
  p <- length(ar)
  q <- length(ma)
  rev.ar <- rev(ar)
  rev.ma <- rev(ma)
  
  if(!is.null(prev.x) & !is.null(prev.e)) {
    burn.period <- 0
    all.e <- c(utils::tail(prev.e, q), stats::rnorm(n))
    all.x <- c(utils::tail(prev.x, p), rep(0, n))
  } else{
    burn.period <- max(burn.period, 100)
    all.e <- c(stats::rnorm(n+burn.period+p))
    all.x <- c(rep(0,n+burn.period+q))
  }
  for(i in 1:(n+burn.period)) {
    all.x[i+p] = sum(all.x[i:(i+p-1)]*rev.ar) + sum(all.e[i:(i+q)]*c(rev.ma, 1))
  }
  return (list(x = all.x[1:n+p+burn.period],
               e = all.e[1:n+q+burn.period]))
}

get.arma.noise <- function(mod, n, cp) {
  segment.len <- get.segment.len(cp = cp, n = n)
  result = list(x = c(), e = c())
  for(i in 1:length(mod)) {
    one.result <- my.arma.sim(mod = mod[[i]], prev.x = result$x, prev.e = result$e, n = segment.len[i])  
    result$x = c(result$x, one.result$x)
    result$e = c(result$e, one.result$e)
  }
  return (result$x)
}

# ========= helpers ===========
get.not.sim.helper <- function(n, cp, start.slope, change.slope, 
                               intercept, jump.size, 
                               init.sigma = NULL, jump.sigma = NULL, 
                               arma.mod = NULL) { 
  mu <- get.jump.signal(n = n, cp = cp, start.slope = start.slope, change.slope = change.slope,
                        intercept = intercept, jump.size = jump.size)
  e <- NULL
  if(!is.null(init.sigma)) {
    sigma <- NULL
    if(!is.null(jump.sigma)) {
      segment.sigma <- cumsum(c(init.sigma, jump.sigma))  
      segment.len <- get.segment.len(cp = cp, n = n)
      sigma <- unlist(lapply(1:length(segment.sigma), function(i) rep(segment.sigma[i], segment.len[i])))  
    } else{
      sigma <- init.sigma
    }
    e <- stats::rnorm(n)*sigma
  } else if(!is.null(arma.mod)) {
    if(length(arma.mod) == 1) {
      e <- stats::arima.sim(model = arma.mod[[1]], n = n, sd = 1)
    } else if(length(arma.mod) == length(cp)+1) {
      e <- get.arma.noise(mod = arma.mod, n = n, cp = cp)
    } else {
      stop("get.not.sim.helper: wrong arma length")
    }
  } else {
    stop("get.not.sim.helper: provide sigma or arma.mod")
  }
  return (list(signal = mu,
               d = mu+e))
}

get.jump.signal <- function(n, cp, start.slope, change.slope, intercept, jump.size) { 
  segment.len <- get.segment.len(cp = cp, n = n)
  segment.slope <- cumsum(c(start.slope, change.slope))
  slope <- unlist(lapply(1:length(segment.slope), function(i) rep(segment.slope[i], segment.len[i])))
  slope[cp+1] <- slope[cp+1] + jump.size
  mu <- cumsum(c(intercept,slope[-1]))
  return (mu)
}

get.segment.len <- function(cp, n) {
  right.end <- c(cp, n)
  segment.len <- diff(c(0, right.end))
  return (segment.len)
}
