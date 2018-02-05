#' @title Estimate the Quantile-optimal Treatment Regime
#' @description Estimate the Quantile-optimal Treatment Regime by inverse probability of weighting
#' 
#' 
#' 
#' @param data a data frame, containing variables in the \code{moPropen} and \code{RegimeClass} and 
#' a component \code{y} as the response.
#' @param regimeClass a formula specifying the class of treatment regimes to search,
#' e.g. if \code{regimeClass = a~x1+x2}, and then this function will search the class of treatment regimes
#' of the form 
#'      \deqn{d(x)=I\left(\beta_0 +\beta_1  x_1 + \beta_2  x_2 > 0\right).
#'        }{d(x)=I(\beta_0 +\beta_1 * x1  + \beta_2 * x2 > 0).}
#'        Polynomial arguments are also supported.
#'          See also 'Details'.
#' @param tau a value between 0 and 1. This is the quantile of interest.
#' @param cl.setup the number of nodes. >1 indicates choosing parallel computing option in 
#'        \code{rgenoud::genoud}. Default is 1.
#' @param moPropen  The propensity score model for the probability of receiving 
#'        treatment level 1.
#'        When \code{moPropen} equals the string "BinaryRandom",  the proportion of observations
#'        receiving treatment level 1 in the sample will be employed
#'        as a good estimate of the probability for each observation.
#'        Otherwise, this argument should be a formula/string, based on which this function
#'        will fit a logistic regression on the treatment level.  e.g. \code{a1~x1}.
#' @param s.tol This is the tolerance level used by \code{genoud}. 
#'              Default is \eqn{10^{-5}} times the difference between
#'              the largest and the smallest value in the observed responses.
#'              This is particularly important when it comes to evaluating \code{it.num}. 
#' @param pop.size an integer with the default set to be 3000. This is the population number for the first generation
#'                in the genetic algorithm (\code{rgenoud::genoud}).
#' @param max logical. If \code{max=TRUE}, it indicates we wish to maximize the marginal
#'        quantile; if \code{max=FALSE}, we wish to minimize the marginal quantile. The default is \code{TRUE}.
#' @param p_level choose between 0,1,2,3 to indicate different levels of output
#'          from the genetic function. Specifically, 0 (minimal printing),
#'            1 (normal), 2 (detailed), and 3 (debug.)
#' @param it.num integer > 1. This argument will be used in \code{rgeound::geound} function.
#'            If there is no improvement in the objective function in this number of generations,
#'        \code{rgenoud::genoud} will think that it has found the optimum.
#' @param  hard_limit logical. When it is true the maximum number of generations
#'          in  \code{rgeound::geound} cannot exceed 100. Otherwise, in this function, only
#' \code{it.num} softly controls when \code{genoud} stops. Default is \code{FALSE}.
#'
#'
#' @return
#' This function returns an object with 7 objects. Both \code{coefficients}
#' and \code{coef.orgn.scale}  were normalized to have unit euclidean norm.
#' 
#'\describe{
#'  \item{\code{coefficients}}{the parameters indexing the estimated 
#'  quantile-optimal treatment regime for 
#'  standardized covariates.}
#'  \item{\code{coef.orgn.scale}}{the parameter indexing the estimated 
#'  quantile-optimal treatment regime for the original input covariates.}
#'  \item{\code{tau}}{the quantile of interest}
#'  \item{\code{hatQ}}{the estimated marginal tau-th quantile when the treatment 
#'          regime indexed by \code{coef.orgn.scale} is applied on everyone.
#'           See the 'details' for connection between \code{coef.orgn.scale} and
#'           \code{coefficient}.}
#'  \item{\code{call}}{the user's call.}
#'  \item{\code{moPropen}}{the user specified propensity score model}
#'  \item{\code{regimeClass}}{the user specified class of treatment regimes}
#'  }
#' 
#' @details  Note that all estimation functions in this package use the same type
#' of standardization on covariates. Doing so would allow us to provide a bounded 
#' domain of parameters for searching in the genetic algorithm.
#'  
#' This estimated parameters indexing the quantile-optimal treatment regime are returned \emph{in two scales:}
#' \enumerate{
#'    \item The returned \code{coefficients} is the set of parameters after covariates \eqn{X} 
#'    are standardized to be in the interval [0, 1]. To be exact, every covariate is 
#'    subtracted by the smallest observed value and divided by the difference between 
#'    the largest and the smallest value.  Next, we carried out the algorithm in Wang et al. 2017 to get the estimated
#'    regime parameters, \code{coefficients}, based on the standardized data. 
#'    For the identifiability issue, we force the Euclidean norm of \code{coefficients}
#'    to be 1.
#'
#'    \item In contrast, \code{coef.orgn.scale} corresponds to the original covariates,
#'     so the associated decision rule can be applied directly to novel observations. 
#'     In other words, let \eqn{\beta} denote the estimated parameter in the original 
#'    scale, then the estimated treatment regime is:  
#'        \deqn{ d(x)= I\{\hat{\beta}_0 + \hat{\beta}_1 x_1 + ... + \hat{\beta}_k x_k > 0\}.}{
#'         d(x)= I{\beta_0 + \beta_1*x_1 + ... + \beta_k*x_k > 0}.}
#'    The estimated \eqn{\bm{\hat{\beta}}}{\beta} is returned as \code{coef.orgn.scale}.
#'    The same as \code{coefficients}, we force the Euclidean norm of \code{coef.orgn.scale}
#'    to be 1.
#' }
#'     If, for every input covariate, the smallest observed value is exactly 0 and the range 
#'    (i.e. the largest number minus the smallest number) is exactly 1, then the estimated 
#'    \code{coefficients} and \code{coef.orgn.scale} will render identical.
#' 
#' @references 
#' \insertRef{wang2017quantile}{quantoptr}
#' 
#' @author Yu Zhou, \email{zhou0269@umn.edu} with substantial contribution from Ben Sherwood.
#' @export
#' @importFrom rgenoud genoud
#' @import stats
#' @import quantreg
#' @importFrom stringr str_replace_all
#' @importFrom methods is
#'
#'
#' @examples
#' GenerateData <- function(n)
#' {
#'   x1 <- runif(n, min=-0.5,max=0.5)
#'   x2 <- runif(n, min=-0.5,max=0.5)
#'   error <- rnorm(n, sd= 0.5)
#'   tp <- exp(-1+1*(x1+x2))/(1+exp(-1+1*(x1+x2)))
#'   a <- rbinom(n = n, size = 1, prob=tp)
#'   y <-  1+x1+x2 +  a*(3 - 2.5*x1 - 2.5*x2) +  (0.5 + a*(1+x1+x2)) * error
#'   return(data.frame(x1=x1,x2=x2,a=a,y=y))
#' }
#' n <- 300
#' testData <- GenerateData(n)
#'
#' # 1. Estimate the 0.25th-quantile optimal treatment regime. ###
#' \donttest{
#' fit1 <- IPWE_Qopt(data = testData, regimeClass = "a~x1+x2",
#'            tau = 0.25, moPropen="a~x1+x2")
#' fit1
#' }
#'
#' # 2. Go parallel. This saves time in calculation. ###
#' \donttest{
#' fit2 <- IPWE_Qopt(data = testData, regimeClass = "a~x1+x2",
#'            tau = 0.25, moPropen="a~x1+x2", cl.setup=2)
#' fit2
#' }
#' 
#' \dontshow{
#' set.seed(1100)
#' testData2 <- GenerateData(30)
#' fit2.test <- IPWE_Qopt(data = testData, regimeClass = "a~x1+x2",
#'            tau = 0.25, moPropen="a~x1+x2", cl.setup=1, pop.size=500, it.num=1, 
#'            s.tol=0.3)
#' fit2.test
#' }
#' 
#' # 3. Set a quardratic term in the class #######################
#' \donttest{
#' fit3 <- IPWE_Qopt(data = testData, regimeClass = "a~x1+x2+I(x1^2)",
#'                   tau = 0.25, moPropen="a~x1+x2", pop.size=1000)
#' fit3
#' }
#' 
#' # 4. Set screen prints level. #######################
#' # Set the p_level to be 0, 
#' # then all screen prints from the genetic algorithm will be suppressed.
#' \donttest{
#' fit4 <- IPWE_Qopt(data = testData, regimeClass = "a~x1+x2",
#'            tau = 0.25, moPropen="a~x1+x2", cl.setup=2, p_level=0)
#' fit4
#' }
IPWE_Qopt <- function(data, regimeClass, tau, moPropen="BinaryRandom",
                      max=TRUE, 
                      s.tol, it.num=8, hard_limit=FALSE,
                      cl.setup=1, p_level=1, pop.size=3000){
  call <- match.call()
  if (!is(data, "data.frame")) 
    stop("'data' must be a data frame.")

  if(!("y" %in% names(data)))
    stop("The response variable 'y' must be present in 'data'.")
  
  if(missing(s.tol))
    s.tol <- diff(range(data$y))*1e-05
  
  if(!(tau<1 & tau>0))
    stop("The quanitle of interst, 'tau' must be strictly bigger 
         than 0 and smaller than 1.")

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
                                           FUN = function(x) (x-min(x))/(max(x)-min(x)))[,-1])
  
  if (moPropen =="BinaryRandom"){
    ph <- rep(mean(data[,txname]), nrow(p.data))
  } else {
    moPropen <- as.formula(moPropen)
    logistic.model.tx <- glm(formula = moPropen, data = data, family=binomial)
    ph <- as.vector(logistic.model.tx$fit)
  }
  
  est <- qestimate(tau=tau, x=p.data.scale, y=data$y, max=max,
                   a=data[,txname], prob=ph,
                   p_level=p_level, nvars=nvars,
                   cl.setup = cl.setup, it.num = it.num,
                   s.tol=s.tol, pop.size = pop.size, hard_limit=hard_limit)
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
            tau=tau,
            hatQ = est$hatQ,
            call=call,
            moPropen=moPropen,
            regimeClass=regimeClass)
  return(fit)
}


