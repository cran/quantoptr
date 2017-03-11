#' @title Estimate the Marginal Quantile Given a Specific Treatment Regime
#' @description Estimate the marginal quantile if the entire population follows a 
#' treatment regime indexed by the given parameters.
#' This function supports the \code{\link{qestimate}} function.
#' @inheritParams abso_diff_est
#' @param tau a numeric value between 0 and 1. The quantile level of interest.
#'
#' @export
quant_est <- function(beta,x,y,a,prob,tau){
  #quantile estimator
  g <- as.numeric(x %*% beta > 0)
  c <- a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  wts <- c*wts
  wts[wts<0.05] = 0.05
  suppressWarnings(model <- quantreg::rq(y ~ 1, weights=wts, tau=tau))
  return(coefficients(model)[1])
}
