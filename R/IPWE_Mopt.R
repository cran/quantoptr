#' @title Estimate the Mean-optimal Treatment Regime
#' @description \code{IPWE_Mopt} aims at estimating the treatment regime which
#' maximizes the marginal mean of the potential outcomes.
#' 
#' 
#' 
#' @param max logical. If \code{max=TRUE}, it indicates we wish to maximize the marginal
#' mean; If \code{max=FALSE}, we wish to minimize the marginal mean. The default is \code{TRUE}.
#' @inheritParams IPWE_Qopt
#'
#'
#'
#'
#' @return
#' This function returns an object with 6 objects. Both \code{coefficients}
#' and \code{coef.orgn.scale}  were normalized to have unit euclidean norm.
#'\describe{
#'  \item{\code{coefficients}}{the parameters indexing the estimated 
#'  mean-optimal treatment regime for 
#'  standardized covariates.}
#'  \item{\code{coef.orgn.scale}}{the parameter indexing the estimated 
#'  mean-optimal treatment regime for the original input covariates.}
#'  \item{\code{hatM}}{the estimated marginal mean  when a treatment regime indexed by
#'         \code{coef.orgn.scale} is applied on everyone. See the 'details' for
#'          connection between \code{coef.orgn.scale} and
#'           \code{coefficient}.}
#'  \item{\code{call}}{the user's call.}
#'  \item{\code{moPropen}}{the user specified propensity score model}
#'  \item{\code{regimeClass}}{the user specified class of treatment regimes}
#'  }
#'  
#'  
#'  
#' @details  Note that all estimation functions in this package use the same type
#' of standardization on covariates. Doing so would allow us to provide a bounded 
#' domain of parameters for searching in the genetic algorithm.
#'  
#' This functions returns the estimated parameters indexing the 
#' mean-optimal treatment regime under two scales. 
#' 
#' The returned \code{coefficients} is the set of parameters when covariates are 
#' all standardized to be in the interval [0, 1] by subtracting the smallest observed
#'  value and divided by the difference between the largest and the smallest value. 
#'
#' While the returned \code{coef.orgn.scale} corresponds to the original covariates,
#' so the associated decision rule can be applied directly to novel observations. 
#' In other words, let \eqn{\beta} denote the estimated parameter in the original 
#' scale, then the estimated treatment regime is:  
#' \deqn{ d(x)= I\{\hat{\beta}_0 + \hat{\beta}_1 x_1 + ... + \hat{\beta}_k x_k > 0\}.}{
#'  d(x)= I{\beta_0 + \beta_1*x_1 + ... + \beta_k*x_k > 0}.}
#' The estimated \eqn{\bm{\hat{\beta}}}{\beta} is returned as \code{coef.orgn.scale}.
#' 
#' If, for every input covariate, the smallest observed value is exactly 0 and the range 
#' (i.e. the largest number minus the smallest number) is exactly 1, then the estimated 
#' \code{coefficients} and \code{coef.orgn.scale} will render identical.
#' 
#'
#' @author Yu Zhou, \email{zhou0269@umn.edu}, with substantial contribution from Ben Sherwood.
#' 
#' 
#' @references 
#' \insertRef{zhang2012robust}{quantoptr}
#' 
#' 
#' @export
#' @importFrom rgenoud genoud
#' @import stats
#' @importFrom stringr str_replace_all
#' @importFrom methods is
#'
#' 
#' 
#' 
#' @examples
#' GenerateData.test.IPWE_Mopt <- function(n)
#' {
#'   x1 <- runif(n)
#'   x2 <- runif(n)
#'   tp <- exp(-1+1*(x1+x2))/(1+exp(-1+1*(x1+x2)))
#'   error <- rnorm(length(x1), sd=0.5)
#'   a <- rbinom(n = n, size = 1, prob=tp)
#'   y <- 1+x1+x2 +  a*(3 - 2.5*x1 - 2.5*x2) + 
#'         (0.5 + a*(1+x1+x2)) * error
#'   return(data.frame(x1=x1,x2=x2,a=a,y=y))
#' }
#' \donttest{
#' n <- 500
#' testData <- GenerateData.test.IPWE_Mopt(n)
#' fit <- IPWE_Mopt(data=testData, regimeClass = a~x1+x2, 
#'                  moPropen=a~x1+x2, 
#'                  pop.size=1000)
#' fit
#' }
#' \dontshow{
#' set.seed(1101)
#' testData <- GenerateData.test.IPWE_Mopt(50)
#' fit <- IPWE_Mopt(data = testData, regimeClass = a~x1+x2, moPropen=a~x1+x2, 
#'                  pop.size=500, it.num=2)
#' fit
#' }
#' 
IPWE_Mopt <- function(data, regimeClass,
                      moPropen="BinaryRandom",
                      max=TRUE,
                      s.tol=0.0001,
                      cl.setup=1, p_level=1,
                      it.num=10, hard_limit=FALSE,
                      pop.size=3000){
  call <- match.call()
  if (!is(data, "data.frame")) 
    stop("'data' must be a data frame.")
  
  if(!("y" %in% names(data)))
    stop("The response variable 'y' must be present in 'data'.")
  
  numNAy <- sum(is.na(data$y))
  if (numNAy>0){
    yNA.idx <- which(is.na(data$y))
    data<-data[!is.na(data$y),]
    message(paste("(", numNAy,
                  "observations are removed since outcome is missing)"))
  }

  regimeClass <- as.formula(regimeClass)
  txname <- as.character(regimeClass[[2]])
  txVec <- try(data[, txname], silent = TRUE)
  if (is(txVec, "try-error")) {
    stop("Variable '", paste0(txname, "' not found in 'data'."))
  }
  if(!all(unique(txVec) %in% c(0,1)))
    stop("The levels of treatment must be numeric, being either 0 or 1.")
  
  # extract the names of the covariates in the decision rule
  p.data <- model.matrix(regimeClass, data)
  minVec <- apply(p.data, MARGIN = 2, min)
  spanVec <- apply(p.data, MARGIN = 2, FUN=function(x) max(x)-min(x))
  
  
  # Dimension of the regimeClass
  nvars <- ncol(p.data) 
  # Rescale each nonconstant variable in regimeClass to range between 0 and 1
  p.data.scale <- cbind(Intercept=1, apply(p.data, MARGIN = 2, 
                                           FUN = function(x) ( x-min(x))/(max(x)-min(x)))[,-1])

  if (moPropen =="BinaryRandom"){
    ph <- rep(data[,txname], nrow(p.data))
  } else {
    moPropen <- as.formula(moPropen)
    logistic.model.tx <- glm(formula = moPropen, data = data, family=binomial)
    ph <- as.vector(logistic.model.tx$fit)
  }
  

  est <- mestimate(x=p.data.scale, y=data$y,
                   a=data[,txname], ph,
                   max=max,
                   p_level, nvars, cl.setup = cl.setup,
                   s.tol=s.tol, it.num=it.num,
                   hard_limit=hard_limit,
                   pop.size = pop.size)
  # parameter indexing the estimated optimal treatment regime, where all
  # covariates are scaled to be from 0 to 1
  coefficient <- est$coefficient
  
  # parameter indexing the same estimated optimal treatment regime, where all
  # covariates are in the original scale
  coef.orgn.scale <- rep(0,length(coefficient))
  coef.orgn.scale[1]  <-coefficient[1]-
    sum(coefficient[-1]*minVec[-1]/spanVec[-1])
  coef.orgn.scale[-1] <- coefficient[-1]/spanVec[-1]
  coef.orgn.scale <- scalar1(coef.orgn.scale)
  
  names(coefficient) <- names(coef.orgn.scale)<- colnames(p.data.scale) 

  fit<-list(coefficients = coefficient,
            coef.orgn.scale = coef.orgn.scale,
            hatM = est$hatM,
            call=call,
            moPropen=moPropen,
            regimeClass=regimeClass)
  return(fit)
}





