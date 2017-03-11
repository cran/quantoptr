#' @title Estimation of the Optimal Treatment Regime defined as Minimizing
#' Gini's Mean Differences
#' @description  
#'  \code{IPWE_MADopt} seeks to estimated the treatment regime which \bold{minimizes}
#'  the Gini's Mean difference defined below.
#' 
#' Besides mean and quantile criterion, in some applications
#' people seek minimization of dispersion in the outcome, which, for example, can
#' be described by Gini's mean difference. Formally, it is defined as the absolute
#' differences of two random variables \eqn{Y_1} and \eqn{Y_2} drawn independently
#' from the same distribution: \deqn{MAD:= E(|Y_1-Y_2|).}
#'
#' Given a treatment regime \eqn{d}, define the potential outcome of a subject
#' following the treatment recommended by \code{d} as
#' \eqn{Y^{*}(d)}. When \eqn{d} is followed by everyone in the target population,
#'  the Gini's mean absolute difference is
#'  \deqn{MAD(d):= E(| Y_1^{*}(d)-Y_2^{*}(d) |).}
#'

#' @inheritParams IPWE_Qopt
#' 
#' 
#' @return
#' This function returns an object with 6 objects. Both \code{coefficients}
#' and \code{coef.orgn.scale}  were normalized to have unit euclidean norm.
#'\describe{
#'  \item{\code{coefficients}}{the parameters indexing the estimated 
#'  MAD-optimal treatment regime for 
#'  standardized covariates.}
#'  \item{\code{coef.orgn.scale}}{the parameter indexing the estimated 
#'  MAD-optimal treatment regime for the original input covariates.}
#'  \item{\code{hat_MAD}}{the estimated MAD when a treatment regime indexed by
#'         \code{coef.orgn.scale} is applied on everyone. See the 'details' for
#'          connection between \code{coef.orgn.scale} and
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
#' This estimated parameters indexing the MAD-optimal treatment regime are returned 
#' \emph{in two scales:}
#' \enumerate{
#'    \item The returned \code{coefficients} is the set of parameters after covariates \eqn{X} 
#'    are standardized to be in the interval [0, 1]. To be exact, every covariate is 
#'    subtracted by the smallest observed value and divided by the difference between 
#'    the largest and the smallest value.  Next, we carried out the algorithm in Wang 2016 to get the estimated
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
#' 
#' 
#' @references 
#' \insertRef{wang2016quant}{quantoptr}
#' 
#' @export
#' 
#' @import utils
#' @importFrom rgenoud genoud
#' @importFrom methods is
#' @examples
#' GenerateData.MAD <- function(n)
#' {
#'   x1 <- runif(n)
#'   x2 <- runif(n)
#'   tp <- exp(-1+1*(x1+x2))/(1+exp(-1+1*(x1+x2)))
#'   a<-rbinom(n = n, size = 1, prob=tp)
#'   error <- rnorm(length(x1))
#'   y <- (1 + a*0.6*(-1+x1+x2<0) +  a*-0.6*(-1+x1+x2>0)) * error
#'   return(data.frame(x1=x1,x2=x2,a=a,y=y))
#' }
#' # The true MAD optimal treatment regime for this generative model
#' # can be deduced trivially, and it is:  c( -0.5773503,  0.5773503,  0.5773503).
#' \dontshow{
#'   set.seed(1103)
#'   testData <- GenerateData.MAD(30)
#'   fit0 <- IPWE_MADopt(data = testData, regimeClass = a~x1+x2,
#'              moPropen=a~x1+x2, s.tol=0.2,
#'              pop.size=300, it.num=2)
#' }
#'
#' # With correctly specified propensity model   ####
#' \donttest{
#' n <- 400
#' testData <- GenerateData.MAD(n)
#' fit1 <- IPWE_MADopt(data = testData, regimeClass = a~x1+x2,
#'                     moPropen=a~x1+x2, cl.setup=2)
#' fit1
#' }
#'              
#' 
#'
#' # With incorrectly specified propensity model ####
#' \donttest{
#' fit2 <- IPWE_MADopt(data = testData, regimeClass = a~x1+x2,
#'                     moPropen="BinaryRandom", cl.setup=2)
#' fit2
#' }
#'


IPWE_MADopt<-function(data, regimeClass,
                      moPropen="BinaryRandom",
                      s.tol, it.num=8, hard_limit=FALSE,
                      cl.setup=1, p_level=1, pop.size=3000 ){
  call <- match.call()
  if (!is(data, "data.frame")) 
    stop("'data' must be a data frame.")
  
  if(!("y" %in% names(data)))
    stop("The response variable 'y' must be present in 'data'.")
  
  if(missing(s.tol))
    s.tol <- diff(range(data$y))*1e-05
  
  
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


  # cl.setup is the number of cores to use
  if(cl.setup>1){
    # parallel computing option
    if((get_os()== "windows"))
      clnodes <- parallel::makeCluster(cl.setup, type="PSOCK")
    else if((get_os()== "osx") | (get_os()== "linux"))
      clnodes <- parallel::makeForkCluster(nnodes =getOption("mc.cores",cl.setup))
    else
      # no parallel
      clnodes <- FALSE
  } else {
    # no parallel
    clnodes <- FALSE
  }

  
  # estimation of the MAD optimal treatment regime
  Cnobs <- combn(1:nrow(data), 2)
  Domains <-cbind(rep(-1,nvars),rep(1,nvars))
  est <-genoud(fn=abso_diff_est, nvars=nvars,
               x=p.data.scale, y=data$y, a=data[,txname],
               prob=ph, Cnobs= Cnobs,
               print.level=p_level, max=FALSE,
               pop.size=pop.size, wait.generations=it.num,
               gradient.check=FALSE, BFGS=FALSE,
               P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
               Domains=Domains, hard.generation.limit=hard_limit,
               starting.values=NULL, solution.tolerance=s.tol,
               optim.method="Nelder-Mead", cluster = clnodes)
  if("cluster" %in%class(clnodes)) { parallel::stopCluster(clnodes) }

  # parameter indexing the estimated optimal treatment regime, where all
  # covariates are scaled to be from 0 to 1
  coefficient <- scalar1(est$par)
  
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
            hat_MAD = est$value,
            call=call,
            moPropen=moPropen,
            regimeClass=regimeClass)
  

  return(fit)
}
