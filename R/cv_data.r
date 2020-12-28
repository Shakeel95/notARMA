get.train.validate.test.data <- function(d, test.n = 10, validate.n =round((length(d)-test.n)*0.15)) {
  n <- length(d)
  train.d <- d[1:(n-validate.n-test.n)]
  validate.d <- d[(n-validate.n-test.n+1):(n-test.n)]
  test.d <- d[(n-test.n+1):n]
  return (list(train.d = train.d,
               validate.d = validate.d,
               test.d = test.d))
}

get.validate.idx <- function(d, include.test.data) {
  if(include.test.data) {
    return (length(d$train.d) + (1:(length(d$validate.d)+length(d$test.d))))
  } else {
    return (length(d$train.d) + (1:length(d$validate.d)))
  }
}

get.test.idx <- function(d) {
  return (length(d$train.d) + length(d$validate.d) + (1:length(d$test.d)))
}

get.test.idx.from.train.n.test.d <- function(train.d, test.d) {
  return (length(train.d) + (1:length(test.d)))
}