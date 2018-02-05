#' @title The Doubly Robust Estimator of the Quantile-Optimal Treatment Regime
#' @description \code{DR_Qopt} implements the doubly robust estimation method to
#' estimate the quantile-optimal treatment regime. The double robustness
#'  property means that it is consistent when either the propensity score model 
#'  is correctly specified, or the conditional quantile function is correctly specified.
#'  Both linear and nonlinear conditional quantile models are considered. See 'Examples'
#'  for an illustrative example.
#' 
#' 
#' 
#' @inheritParams IPWE_Qopt
#' @param data a data frame, must contain all the variables that appear in \code{moPropen},
#'        \code{RegimeClass},  \code{moCondQuant_0},  \code{moCondQuant_1}, and a column named
#'        \code{y} as the observed response. 
#' @param nlCondQuant_0 Logical. When \code{nlCondQuant_0=TRUE},
#' this means the prespecified model for
#' the conditional quantile function given a=0 is nonlinear,
#' so the provided \code{moCondQuant_0}
#' should be nonlinear.
#' @param nlCondQuant_1 Logical. When \code{nlCondQuant_1=TRUE},
#' this means the prespecified model for the conditional quantile function
#' given a=1 is nonlinear,
#' so the provided \code{moCondQuant_1}
#' should be nonlinear.
#' @param length.out an integer greater than 1.  If one of the conditional quantile
#'  model is set to be nonlinear, this argument will be triggered and we will fit 
#'  \code{length.out} models across quantiles equally spaced between 0.001 and 0.999.
#'  Default is 200.
#' @param moCondQuant_0 Either a formula or a string representing
#' the parametric form of the conditional quantile function given that treatment=0.
#' @param moCondQuant_1 Either a formula or a string representing
#' the parametric form of the conditional quantile function given that treatment=1.
#' @param start_0 a named list or named numeric vector of starting estimates for
#' the conditional quantile function when \code{treatment = 0}. This is required when
#' \code{nlCondQuant_0=TRUE}.
#' @param start_1 a named list or named numeric vector of starting estimates for
#' the conditional quantile function when \code{treatment = 1}. This is required when
#' \code{nlCondQuant_1=TRUE}.
#' 
#' 
#' @return 
#' This function returns an object with 9 objects. Both \code{coefficients}
#' and \code{coef.orgn.scale}  were normalized to have unit euclidean norm.
#' \describe{
#'  \item{\code{coefficients}}{the parameters indexing the estimated 
#'  quantile-optimal treatment regime for 
#'  standardized covariates. }
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
#'  \item{\code{moCondQuant_0}}{the user specified conditional quantile model for treatment 0}
#'  \item{\code{moCondQuant_1}}{the user specified conditional quantile model for treatment 1}
#'  }
#' 
#'
#'
#'
#' @details  
#' \itemize{
#'    \item Standardization on covariates AND explanation on the differences between
#'    the two returned regime parameters.
#'        
#'        Note that all estimation functions in this package use the same type
#'        of standardization on covariates. Doing so would allow us to provide a bounded 
#'         domain of parameters for searching in the genetic algorithm.
#'  
#'        This estimated parameters indexing the quantile-optimal treatment regime are returned \emph{in two scales:}
#'   \enumerate{
#'       \item The returned \code{coefficients} is the set of parameters after covariates \eqn{X} 
#'    are standardized to be in the interval [0, 1]. To be exact, every covariate is 
#'    subtracted by the smallest observed value and divided by the difference between 
#'    the largest and the smallest value.  Next, we carried out the algorithm in Wang 2016 to get the estimated
#'    regime parameters, \code{coefficients}, based on the standardized data. 
#'    For the identifiability issue, we force the Euclidean norm of \code{coefficients}
#'    to be 1.
#'
#'        \item In contrast, \code{coef.orgn.scale} corresponds to the original covariates,
#'     so the associated decision rule can be applied directly to novel observations. 
#'     In other words, let \eqn{\beta} denote the estimated parameter in the original 
#'    scale, then the estimated treatment regime is:  
#'        \deqn{ d(x)= I\{\hat{\beta}_0 + \hat{\beta}_1 x_1 + ... + \hat{\beta}_k x_k > 0\}.}{
#'         d(x)= I{\beta_0 + \beta_1*x_1 + ... + \beta_k*x_k > 0}.}
#'    The estimated \eqn{\bm{\hat{\beta}}}{\beta} is returned as \code{coef.orgn.scale}.
#'    The same as \code{coefficients}, we force the Euclidean norm of \code{coef.orgn.scale}
#'    to be 1.
#'    }
#'     If, for each input covariate, the smallest observed value is exactly 0 and the range 
#'    (i.e. the largest number minus the smallest number) is exactly 1, then the estimated 
#'    \code{coefficients} and \code{coef.orgn.scale} will render identical.
#'    
#'    
#'    
#'    \item Property of the doubly robust(DR) estimator. The DR estimator \code{DR_Qopt}
#'    is consistent if either the propensity score model or the conditional quantile
#'    regression model is correctly specified. (Wang et. al. 2016)
#' }
#' 
#' 
#' 
#' @seealso  \code{\link{dr_quant_est}}, \code{\link{augX}}
#' 
#' @author Yu Zhou, \email{zhou0269@umn.edu}
#' @export
#'
#' @references 
#' \insertRef{wang2017quantile}{quantoptr}
#' 
#'
#'
#' @importFrom  rgenoud genoud
#' @import stats
#' @import quantreg
#' @importFrom stringr str_replace_all
#' @examples
#' ilogit <- function(x) exp(x)/(1 + exp(x))
#' GenerateData.DR <- function(n)
#' {
#'  x1 <- runif(n,min=-1.5,max=1.5)
#'  x2 <- runif(n,min=-1.5,max=1.5)
#'  tp <- ilogit( 1 - 1*x1^2 - 1* x2^2)
#'  a <-rbinom(n,1,tp)
#'  y <- a * exp(0.11 - x1- x2) + x1^2 + x2^2 +  a*rgamma(n, shape=2*x1+3, scale = 1) +
#'  (1-a)*rnorm(n, mean = 2*x1 + 3, sd = 0.5)
#'  return(data.frame(x1=x1,x2=x2,a=a,y=y))
#' }

#'
#' regimeClass <- as.formula(a ~ x1+x2)
#' moCondQuant_0 <- as.formula(y ~ x1+x2+I(x1^2)+I(x2^2))
#' moCondQuant_1 <- as.formula(y ~ exp( 0.11 - x1 - x2)+ x1^2 + p0 + p1*x1
#'                            + p2*x1^2 + p3*x1^3 +p4*x1^4 )
#' start_1 = list(p0=0, p1=1.5, p2=1, p3 =0,p4=0)
#' 

#'\dontshow{   
#'   n.test<-30
#'   set.seed(1200)
#'   testdata2 <- GenerateData.DR(n.test)
#'   fit0 <- DR_Qopt(data=testdata2, regimeClass = a ~ x1+x2, tau = 0.2,
#'                  moPropen = a~I(x1^2)+I(x2^2),
#'                  moCondQuant_0 = moCondQuant_0,
#'                  moCondQuant_1 = moCondQuant_1,
#'                  length.out = 2, 
#'                  p_level=1, s.tol=0.5,
#'                  nlCondQuant_1 = TRUE,  start_1=start_1,
#'                  pop.size = 500, it.num =1)
#'}
#'
#' n <- 400
#' testdata <- GenerateData.DR(n)
#'
#' ## Examples below correctly specified both the propensity model and 
#' ##  the conditional quantile model.
#'  \donttest{ 
#'  system.time(
#'  fit1 <- DR_Qopt(data=testdata, regimeClass = regimeClass, 
#'                  tau = 0.25,
#'                  moPropen = a~I(x1^2)+I(x2^2),
#'                  moCondQuant_0 = moCondQuant_0,
#'                  moCondQuant_1 = moCondQuant_1,
#'                  nlCondQuant_1 = TRUE,  start_1=start_1,
#'                  pop.size = 1000))
#'  fit1}
#'  ## Go parallel for the same fit. It would save a lot of time.
#'  ### Could even change the cl.setup to larger values 
#'  ### if more cores are available.
#'  \donttest{ 
#'  system.time(fit2 <- DR_Qopt(data=testdata, regimeClass = regimeClass, 
#'                  tau = 0.25,
#'                  moPropen = a~I(x1^2)+I(x2^2),
#'                  moCondQuant_0 = moCondQuant_0,
#'                  moCondQuant_1 = moCondQuant_1,
#'                  nlCondQuant_1 = TRUE,  start_1=start_1,
#'                  pop.size = 1000, cl.setup=2))
#'  fit2}
#'




DR_Qopt<-function(data, regimeClass,
                  tau, moPropen = "BinaryRandom",
                  nlCondQuant_0=FALSE, nlCondQuant_1=FALSE,
                  moCondQuant_0,
                  moCondQuant_1,
                  max=TRUE,
                  length.out=200,
                  s.tol, 
                  it.num = 8,
                  cl.setup=1, p_level=1, pop.size=3000,
                  hard_limit=FALSE, start_0=NULL, start_1=NULL)
{
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
  foo <- augX(raw.data=data, 
              length.out = length.out, txVec = txVec,
              moCondQuant_0=moCondQuant_0,
              moCondQuant_1=moCondQuant_1,
              nlCondQuant_0=nlCondQuant_0,
              nlCondQuant_1=nlCondQuant_1,
              start_0=start_0, start_1=start_1, 
              clnodes=clnodes)
  y.a.0<-foo$y.a.0
  y.a.1<-foo$y.a.1

  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  # the standardized covariates are used in optimization
  est <- genoud(fn=dr_quant_est, nvars=nvars,
              x=p.data.scale, 
              y=data$y, a=data[,txname], prob=ph, tau=tau,
              y.a.0=y.a.0, y.a.1=y.a.1,
              max=TRUE, print.level=p_level, pop.size=pop.size,
              wait.generations=it.num,gradient.check=FALSE, BFGS=FALSE,
              P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              Domains=Domains, starting.values=rep(0, nvars),
              hard.generation.limit=hard_limit,
              solution.tolerance=s.tol, optim.method="Nelder-Mead",
              cluster = clnodes)
  if("cluster" %in% class(clnodes)) { parallel::stopCluster(clnodes) }

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
  
  ## see number of minimum values
  nMinimum <- dr_quant_est( beta = coefficient, 
                            x= p.data.scale, 
                            y=data$y, a=txVec,
                            prob=ph, tau=tau,
                            y.a.0=y.a.0,
                            y.a.1=y.a.1, 
                            num_min =TRUE)
  
  names(coefficient) <- names(coef.orgn.scale) <- colnames(p.data.scale) 
  fit<-list(coefficients = coefficient,
            coef.orgn.scale = coef.orgn.scale,
            tau=tau,
            hatQ=est$value,
            call=call,
            moPropen=moPropen,
            regimeClass=regimeClass,
            nMinimum=nMinimum,
            moCondQuant_0 = moCondQuant_0,
            moCondQuant_1 = moCondQuant_1)

  return(fit)
}

