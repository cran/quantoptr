#' @title
#' The Inverse Probability Weighted Estimator of the Marginal Mean Given a Specific Treatment Regime
#' @description
#' Estimate the marginal mean of the response when the entire population
#' follows a treatment regime. This function implements the inverse probability weighted
#' estimator proposed by Baqun Zhang et. al..
#' 
#' This function supports the \code{\link{mestimate}} function.
#'
#' @inheritParams abso_diff_est
#'
#' @references 
#' \insertRef{zhang2012robust}{quantoptr}
#' 
#' @export
mean_est<-function(beta,x,a,y,prob)
{
  #mean estimator
  g<-as.numeric(x%*%beta>0)
  c<-a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  val <- mean(c*y*wts)

  return(val)
}
