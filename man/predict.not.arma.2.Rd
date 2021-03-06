% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/not_arma_2.r
\name{predict.not.arma.2}
\alias{predict.not.arma.2}
\title{one step ahead prediction form fitted NOT-ARMA-2 model}
\usage{
\method{predict}{not.arma.2}(object, test.i, ...)
}
\arguments{
\item{object}{fitted model from not.arma.2()}

\item{test.i}{time point to predict}

\item{...}{additional arguments}
}
\value{
one step ahead value
}
\description{
one step ahead prediction form fitted NOT-ARMA-2 model
}
\examples{
\donttest{
d <- get.s1.sim(n = 200, cp = 100, intercept = 0, jump.size = 1, 
sigma = 1, arma.mod = c(0.8, 0.2))$d
mod <- not_arma_2(d[1:199], contrast = "pcwsConstMean")
pred <- predict(mod, 200)
}
}
