#' @title Estimate the Two-stage Quantile-optimal Treatment Regime
#' @description  This function implements the estimator of 
#' two-stage quantile-optimal treatment regime by inverse probability of weighting 
#' proposed by Lan Wang, et al. As there are more than one stage, the second stage
#' treatment regime could take into account the evolving status of an individual
#' after the first stage
#' and the treatment level received in the first stage. We assume the options at 
#' the two stages are both binary and take the form:
#' \deqn{d_1(x)=I\left(\beta_{10} +\beta_{11}  x_{11} +...+ \beta_{1k}  x_{1k} > 0\right), 
#' }{
#' d1(x)=I(beta_10 +beta_11 * x11 +... + beta_1k * x1k > 0),
#' }
#' \deqn{d_2(x)=I\left(\beta_{20} +\beta_{21}  x_{21} +...+ \beta_{2p}  x_{2p} > 0\right)}{
#' d2(x)=I(beta_20 +beta_21 * x21 +... + beta_2p * x2p > 0) }
#' 
#'
#' @inheritParams IPWE_Qopt
#' @inheritParams TwoStg_Mopt
#' 
#' 
#' 
#' 
#' @return
#' This function returns an object with 7 objects. Both \code{coefficients}
#' and \code{coef.orgn.scale}  were normalized to have unit euclidean norm.
#'\describe{
#'  \item{\code{coef.1}, \code{coef.2}}{the set of parameters indexing the estimated 
#'  quantile-optimal treatment regime for 
#'  standardized covariates.}
#'  \item{\code{coef.orgn.scale.1}, \code{coef.orgn.scale.2}}{the set of parameter 
#'  indexing the estimated quantile-optimal treatment regime for the original input covariates.}
#'  \item{\code{tau}}{the quantile of interest}
#'  \item{\code{hatQ}}{the estimated marginal quantile when the treatment 
#'          regime indexed by \code{coef.orgn.scale.1} and \code{coef.orgn.scale.2} 
#'          is applied on the entire population.
#'           See the 'details' for connection between \code{coef.orgn.scale.k} and
#'           \code{coef.k}.}
#'  \item{\code{call}}{the user's call.}
#'  \item{\code{moPropen1}, \code{moPropen2}}{the user specified propensity score models
#'  for the first and the second stage respectively}
#'  \item{\code{regimeClass.stg1},  \code{regimeClass.stg2}}{the user specified 
#'  class of treatment regimes for the first and the second stage respectively}
#'}
#'
#'
#' @details
#'   Note that all estimation functions in this package use the same type
#' of standardization on covariates. Doing so would allow us to provide a bounded 
#' domain of parameters for searching in the genetic algorithm.
#'  
#' 
#' For every stage \code{k}, \eqn{k=1,2}, this estimated parameters indexing the 
#' two-stage quantile-optimal treatment regime are returned \emph{in two scales:}
#' \enumerate{
#'    \item , the returned \code{coef.k} 
#'      is the set of parameters that we estimated after standarding
#'      every covariate available for decision-making
#'      at stage \code{k} to be in the interval [0, 1]. To be exact, every covariate is 
#'    subtracted by the smallest observed value and divided by the difference between 
#'    the largest and the smallest value.  Next, we carried out the algorithm in Wang et. al. 2016 to get the estimated
#'    regime parameters, \code{coef.k}, based on the standardized data. 
#'    For the identifiability issue, we force the Euclidean norm of \code{coef.k}
#'    to be 1.
#'
#'  \item The difference between \code{coef.k} and \code{coef.orgn.scale.k} is that the latter
#'    set of parameters correspond to the original covariates,
#'    so the associated decision rule can be applied directly to novel observations. 
#'    In other words, let \eqn{\beta} denote the estimated parameter in the original 
#'    scale, then the estimated treatment regime is:  
#'      \deqn{ d(x)= I\{\beta_0 + \beta_1 x_1 + ... + \beta_k x_k > 0\},}{ d(x)=
#'        I{\beta_0 + \beta_1*x_1 + ... + \beta_k*x_k > 0},}
#'    where the \eqn{\beta} values are returned as \code{coef.orgn.scale.k}, and the
#'    the vector \eqn{(1, x_1,...,x_k)} corresponds to the specified class of treatment
#'    regimes in the \code{k}th stage.
#' }
#' 
#' If, for every input covariate, the smallest observed value is exactly 0 and the range 
#' (i.e. the largest number minus the smallest number) is exactly 1, then the estimated 
#' \code{coef.k} and \code{coef.orgn.scale.k} will render identical.

#' 
#' 
#' @references 
#' \insertRef{wang2016quant}{quantoptr}
#' 
#' @author Yu Zhou, \email{zhou0269@umn.edu}
#' @export
#' @importFrom  rgenoud genoud
#' @import stats
#' @import quantreg
#' @examples
#' library(faraway)
#' GenerateData.2stg <- function(n){
#'  x1 <- runif(n)
#'  p1 <- ilogit(-0.5+x1)
#'  a1 <- rbinom(n, size=1, prob=p1)
#'  
#'  x2 <- runif(n,x1,x1+1)
#'  p2 <- ilogit(-1 + x2)
#'  a2 <- rbinom(n, size=1, prob=p2)
#'  
#'  mean <- 1+x1+a1*(1-3*(x1-0.2)^2) +x2 + a2*(1-x2-x1)
#'  y <- mean + (1+a1*(x1-0.5)+0.5*a2*(x2-1))*rnorm(n,0,sd = 1)
#'  return(data.frame(x1,a1,x2,a2,y))
#' }
#' \donttest{
#' n <- 400
#' testdata <- GenerateData.2stg(n)
#' fit <- TwoStg_Qopt(data=testdata, tau=0.2,
#'                    regimeClass.stg1=a1~x1, regimeClass.stg2=a2~x1+a1+x2,
#'                    moPropen1=a1~x1, moPropen2=a2 ~ x2,
#'                    cl.setup=2)
#' fit
#' }
#' 
#' \dontshow{
#' set.seed(11001)
#' testdata <- GenerateData.2stg(30)
#' fit <- TwoStg_Qopt(data=testdata, tau=0.2, 
#'                    regimeClass.stg1=a1~x1, regimeClass.stg2=a2~x2,
#'                    moPropen1=a1~x1, moPropen2=a2 ~ x2,
#'                    p_level=0,
#'                    cl.setup=1, pop.size=300, it.num=1, s.tol=0.1)
#' }

TwoStg_Qopt <- function(data, tau,
                        regimeClass.stg1,
                        regimeClass.stg2,
                        moPropen1="BinaryRandom",
                        moPropen2="BinaryRandom",
                        s.tol=0.0001, it.num=8,
                        max=TRUE,
                        cl.setup=1, p_level=1, pop.size=1000, hard_limit=FALSE){
  call <- match.call()
  if (!is(data, "data.frame")) 
    stop("'data' must be a data frame.")
  
  if(!(tau<1 & tau>0))
    stop("The quanitle of interst, 'tau' must be strictly bigger 
         than 0 and smaller than 1.")
  
  if(!("y" %in% names(data)))
    stop("The response variable 'y' must be present in 'data'.")
  
  numNAy <- sum(is.na(data$y))
  if (numNAy>0){
    yNA.idx <- which(is.na(data$y))
    data<-data[!is.na(data$y),]
    message(paste("(", numNAy,
                  "observations are removed since outcome is missing)"))
  }
  
  regimeClass.stg1 <- as.formula(regimeClass.stg1)
  regimeClass.stg2 <- as.formula(regimeClass.stg2)
  # extract the names of the covariates in the decision rule
  p.data1 <- model.matrix(regimeClass.stg1, data)
  p.data2 <- model.matrix(regimeClass.stg2, data)
  txname.stg1 <- as.character(regimeClass.stg1[[2]])
  txname.stg2 <- as.character(regimeClass.stg2[[2]])
  

  minVec1 <- apply(p.data1, MARGIN = 2, min)
  spanVec1 <- apply(p.data1, MARGIN = 2, FUN=function(x) max(x)-min(x))
  
  minVec2 <- apply(p.data2, MARGIN = 2, min)
  spanVec2 <- apply(p.data2, MARGIN = 2, FUN=function(x) max(x)-min(x))
  
  # Rescale each nonconstant variable in regimeClass to range between 0 and 1
  p.data.scale1 <- cbind(Intercept=1, apply(p.data1, MARGIN = 2, 
                                            FUN = function(x) (x-min(x))/(max(x)-min(x)))[,-1])
  p.data.scale2 <- cbind(Intercept=1, apply(p.data2, MARGIN = 2, 
                                            FUN = function(x) (x-min(x))/(max(x)-min(x)))[,-1])
  
  
  txVec1 <- try(data[, txname.stg1], silent = TRUE)
  if (is(txVec1, "try-error")) {
    stop("Variable '", paste0(txname.stg1, "' not found in 'data'."))
  }
  if(!all(unique(txVec1) %in% c(0,1)))
    stop("The levels of treatment in the first stage must be numeric, 
         being either 0 or 1.")
  
  txVec2 <- try(data[, txname.stg2], silent = TRUE)
  if (is(txVec2, "try-error")) {
    stop("Variable '", paste0(txname.stg2, "' not found in 'data'."))
  }
  if(!all(unique(txVec2) %in% c(0,1)))
    stop("The levels of treatment in the second stage must be numeric, 
         being either 0 or 1.")
  
  
  nvars.stg1 <- ncol(p.data1)
  nvars.stg2 <- ncol(p.data2)
  nvars.total <- nvars.stg1 + nvars.stg2
  
  
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
  


  if (moPropen1 =="BinaryRandom"){
    ph.stg1 <- rep(mean(txVec1), nrow(p.data1))
  } else {
    moPropen1 <- as.formula(moPropen1)
    logistic.model.tx.stg1 <- glm(moPropen1, data = data, family=binomial)
    ph.stg1 <- as.vector(logistic.model.tx.stg1$fit)
  }

  if (moPropen2 =="BinaryRandom"){
    ph.stg2 <- rep(mean(txVec2), nrow(p.data2))
  } else {
    moPropen2 <- as.formula(moPropen2)
    logistic.model.tx.stg2 <- glm(moPropen2, data = data, family=binomial)
    ph.stg2 <- as.vector(logistic.model.tx.stg2$fit)
  }



  # define the form of decision rules
  CovSpace_1 <- p.data.scale1
  CovSpace_2 <- p.data.scale2 
  
  # calculate the probability of selection
  a1 <- data[,txname.stg1]
  a2 <- data[,txname.stg2]
  prob <- (ph.stg1^a1)*((1-ph.stg1)^(1-a1)) *(ph.stg2^a2)*((1-ph.stg2)^(1-a2))


  Qtauhat_TwoStg <- function(beta12, nv.stg1, CovSpace_1, CovSpace_2,
                            a1, a2, y, prob, tau){
    beta1 <- beta12[1:nv.stg1]
    beta2 <- beta12[-c(1:nv.stg1)]
    g1 <- as.numeric(CovSpace_1%*%beta1 > 0)
    g2 <- as.numeric(CovSpace_2%*%beta2 > 0)
    c_infty <- as.numeric(a1==g1 & a2 == g2)
    wts <- c_infty *( 1/prob )
    model <- rq(y ~ 1, weights=wts, tau=tau)

    criter <- coefficients(model)[1]
    return(criter)
  }
  Domains <- cbind(rep(-1,nvars.total),rep(1,nvars.total))
  est <- genoud(fn=Qtauhat_TwoStg, nvars=nvars.total,
                nv.stg1=nvars.stg1,
                CovSpace_1=CovSpace_1, 
                CovSpace_2=CovSpace_2,
                a1=a1, a2=a2,
                y=data$y, prob=prob, tau=tau,
                print.level=p_level, max=max,
                pop.size=pop.size,
                wait.generations=it.num,
                gradient.check=FALSE, BFGS=FALSE,
                P1=50, P2=50, P3=50, P4=50, P5=50,
                P6=50, P7=50, P8=50, P9=0,
                Domains=Domains,
                starting.values=NULL,
                hard.generation.limit=hard_limit,
                solution.tolerance=s.tol,
                optim.method="Nelder-Mead",
                cluster = clnodes)
  if("cluster" %in%class(clnodes)) { parallel::stopCluster(clnodes) }

  coef.1 <- scalar1(est$par[1:nvars.stg1])
  coef.2 <- scalar1(est$par[-c(1:nvars.stg1)])
  
  # parameter indexing the same estimated optimal treatment regime, where all
  # covariates are in the original scale
  coef.orgn.scale.1 <- rep(0,length(coef.1))
  coef.orgn.scale.2 <- rep(0,length(coef.2))
  
  coef.orgn.scale.1[1] <-coef.1[1]- sum(coef.1[-1]*minVec1[-1]/spanVec1[-1])
  coef.orgn.scale.2[1] <-coef.2[1]- sum(coef.2[-1]*minVec2[-1]/spanVec2[-1])
  
  coef.orgn.scale.1[-1] <- coef.1[-1]/spanVec1[-1]
  coef.orgn.scale.2[-1] <- coef.2[-1]/spanVec2[-1]
  
  coef.orgn.scale.1 <- scalar1(coef.orgn.scale.1)
  coef.orgn.scale.2 <- scalar1(coef.orgn.scale.2)
  
  names(coef.1) <- names(coef.orgn.scale.1)<- colnames(p.data.scale1) 
  names(coef.2) <- names(coef.orgn.scale.2)<- colnames(p.data.scale2)
  fit<- list(coef.1 = coef.1, 
             coef.orgn.scale.1 =coef.orgn.scale.1,
             coef.2=coef.2,
             coef.orgn.scale.2 =coef.orgn.scale.2,
             hatQ = est$value,
             call=call,
             moPropen1=moPropen1, 
             moPropen2=moPropen2,
             regimeClass.stg1 = regimeClass.stg1,
             regimeClass.stg2 = regimeClass.stg2)
  
  return(fit)
}

