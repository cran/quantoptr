### for R_Quant_Est###
# two versions:
# beta_rq, and beta_nlrq are different because the process from rq explicitly states that
# the value of coefficients is the same for tau in [tau_t, tau_t+1)
## Tools for DR method #############################


#' @title Generate Pseudo-Responses Based on Conditional Quantile Regression
#' Models
#' 
#' @description This function supports the \code{\link{DR_Qopt}} function.
#' For every observation, we generate pseudo-observations corresponding
#' to treatment 0 and 1 respectively based on working conditional quantile models.
#' 
#' @param raw.data A data frame, must contain all the variables that appear in 
#' \code{moCondQuant_0} and \code{moCondQuant_1}.
#' @param length.out an integer greater than 1.  If one of the conditional quantile
#'  model is set to be nonlinear, this argument will be triggered and we will fit 
#'  \code{length.out} models across quantiles equally spaced between 0.001 and 0.999.
#'  The larger this value, the more refined the performance of this method.
#'  Default is 200.
#' @param txVec a numeric vector of observed treatment levels coded 0L and 1L.
#' @param moCondQuant_0 A formula, used to specify the formula for the conditional
#'        quantile function when treatment = 0.
#' @param moCondQuant_1 A formula, used to specify the formula for the conditional
#'        quantile function when treatment = 1.
#' @param nlCondQuant_0 logical.
#'        When \code{nlCondQuant_0 = TRUE}, it is indicated that \code{moCondQuant_0} is nonlinear.
#'        The default value of this variable is \code{FALSE}.
#' @param nlCondQuant_1 logical.
#'        When \code{nlCondQuant_1 = TRUE}, it is indicated that \code{moCondQuant_1} is nonlinear.
#'        The default value of this variable is \code{FALSE}.
#'
#' @param start_0 either a list object, providing the starting value in estimating 
#'        the parameters in the nonlinear conditional quantile model, given that treatment=0. 
#'        Default is \code{NULL}, corresponding to the case when \code{nlCondQuant_0=FALSE}.
#' @param start_1 either a list object, providing the starting value in estimating 
#'          the parameters in the nonlinear conditional quantile model, given that treatment=0. 
#'          Default is \code{NULL}, corresponding to the case when \code{nlCondQuant_1=FALSE}.

#' @param clnodes Either a cluster object to enable parallel computation or 
#'        \code{NULL}. If \code{NULL}, no parallel computation will be used.
#'
#' @details 
#' This function implements the algorithm to generate individual level pseudo responses 
#' for two treatment levels respectively.
#' 
#' For each observation, two independent random variables from 
#' \eqn{{unif}[0,1]}{unif[0,1]} are generated. Denote them by \eqn{u_0}{u0} 
#' and \eqn{u_1}{u1}. Approximately, this function then estimates the \eqn{u_0}{u0}th quantile
#' of this observation were treatment level 0 is applied via the conditional \eqn{u_0}{u0}th quantile regression. 
#' This estimated quantile will be the pseudo-response for treatment 0. 
#' Similarly, this function the pseudo-response for treatment 1 will be estimated and returned.
#' 
#' See the reference paper for a more formal explanation.
#' 
#' @references 
#' \insertRef{wang2017quantile}{quantoptr}
#'
#' @return
#' It returns a list object, consisting of the following elements:
#' \enumerate{
#'    \item \code{y.a.0}, the vector of estimated individual level pseudo outcomes,
#'    given the treatment is 0;
#'    \item \code{y.a.1}, the vector of estimated individual level pseudo outcomes, 
#'    given the treatment is 1;
#'    \item \code{nlCondQuant_0}, logical, indicating whether the \code{y.a.0}
#'          is generated based on a nonlinear conditional quantile model.
#'    \item \code{nlCondQuant_1}, logical, indicating whether the \code{y.a.1}
#'          is generated based on a nonlinear conditional quantile model.
#'  }
#'  
#' @examples
#' ilogit <- function(x) exp(x)/(1 + exp(x))
#' GenerateData.DR <- function(n)
#' {
#'   x1 <- runif(n,min=-1.5,max=1.5)
#'   x2 <- runif(n,min=-1.5,max=1.5)
#'   tp <- ilogit( 1 - 1*x1^2 - 1* x2^2)
#'   a <-rbinom(n,1,tp)
#'   y <- a * exp(0.11 - x1- x2) + x1^2 + x2^2 +  a*rgamma(n, shape=2*x1+3, scale = 1) +
#'        (1-a)*rnorm(n, mean = 2*x1 + 3, sd = 0.5)
#'   return(data.frame(x1=x1,x2=x2,a=a,y=y))
#' }
#' regimeClass = as.formula(a ~ x1+x2)
#' moCondQuant_0 = as.formula(y ~ x1+x2+I(x1^2)+I(x2^2))
#' moCondQuant_1 = as.formula(y ~ exp( 0.11 - x1 - x2)+ x1^2 + p0 + p1*x1
#' + p2*x1^2 + p3*x1^3 +p4*x1^4 )
#' start_1 = list(p0=0, p1=1.5, p2=1, p3 =0,p4=0)
#' 
#' \dontrun{
#' n<-200
#' testdata <- GenerateData.DR(n)
#' fit1 <- augX(raw.data=testdata, txVec = testdata$a,
#'              moCondQuant_0=moCondQuant_0, moCondQuant_1=moCondQuant_1,
#'              nlCondQuant_0=FALSE,   nlCondQuant_1=TRUE,
#'              start_1=start_1, 
#'              clnodes=NULL)  
#'  
#'# How to use parallel computing in AugX(): ##
#'  
#'# on Mac OSX/linux
#'  clnodes <- parallel::makeForkCluster(nnodes =getOption("mc.cores",2))
#'  fit2 <- augX(raw.data=testdata, length.out = 5, txVec = testdata$a,
#'              moCondQuant_0=moCondQuant_0, moCondQuant_1=moCondQuant_1,
#'              nlCondQuant_0=FALSE,   nlCondQuant_1=TRUE,
#'              start_1=start_1, 
#'              clnodes=clnodes)  
#'   
#'# on Windows
#'  clnodes <- parallel::makeCluster(2, type="PSOCK")
#'  fit3 <- augX(raw.data=testdata, length.out = 5, txVec = testdata$a,
#'              moCondQuant_0=moCondQuant_0, moCondQuant_1=moCondQuant_1,
#'              nlCondQuant_0=FALSE,   nlCondQuant_1=TRUE,
#'              start_1=start_1, 
#'              clnodes=clnodes)  
#'  }
#' \dontshow{ 
#' n<-100
#' testdata.parallel <- GenerateData.DR(n)
#'   foo <- augX(raw.data=testdata.parallel, length.out = 5, txVec = testdata.parallel$a,
#'              moCondQuant_0=moCondQuant_0, moCondQuant_1=moCondQuant_1,
#'              nlCondQuant_0=FALSE,   nlCondQuant_1=TRUE,
#'              start_1=start_1, clnodes=NULL)
#'              } 
#'  
#'  
#' @export
#' @import quantreg
#' @importFrom parallel parLapply

augX <- function(raw.data, length.out=200,
                 txVec,
                 moCondQuant_0, moCondQuant_1,
                 nlCondQuant_0=FALSE, nlCondQuant_1=FALSE,
                 start_0=NULL, start_1=NULL, clnodes )
{
  data0 <- raw.data[txVec==0,]
  data1 <- raw.data[txVec==1,]
  tauvec <- seq(from=0.001, to=0.999, length.out=length.out)

  if (!nlCondQuant_0){
    # md0 <- rq(moCondQuant_0, data=data0 , tau =2)
    md0.rqsol = rq(moCondQuant_0, data=data0 , tau =2)$sol
    md0tau = md0.rqsol[1,]
    md0tau_beta = rbind(md0tau, md0.rqsol[-(1:3),])
    md0tau_beta = t(md0tau_beta)
  } else {
    # the conditional quantile function is nonlinear
    list_nlrq_0_model <- vec_nlrq_0_tau <- NULL
    # take the parallel computing if possible
    if("cluster" %in% class(clnodes)) {
      # parallel computing
      list_0 =  parLapply(cl=clnodes, X=tauvec,
                          fun=function(tau, moCondQuant_0, data0,start_0){
                                  md0.nls.i = nlrq(formula = moCondQuant_0,
                                             data=data0,
                                             start = start_0, tau = tau)},
                          moCondQuant_0,
                          data0, start_0)

      for(i in 1:length(tauvec)){
        list_nlrq_0_model[[length(list_nlrq_0_model)+1]]<- list_0[[i]]
        vec_nlrq_0_tau <- c(vec_nlrq_0_tau, tauvec[i])
      }
    } else {
      # no parallel
      for(i in 1:length(tauvec)){
        md0.nls.i = nlrq(formula = moCondQuant_0, data=data0,
                         start = start_0, tau = tauvec[i])
        list_nlrq_0_model[[length(list_nlrq_0_model)+1]]<- md0.nls.i
        vec_nlrq_0_tau <- c(vec_nlrq_0_tau, tauvec[i])
      }
    }
  }

  # A=1 ####
  # Next, estimated the conditional quantile regression function for responses
  # given that A=1
  if (!nlCondQuant_1){
    # the conditional quantile function is linear
    md1.rqsol = rq(moCondQuant_1, data=data1 , tau =2)$sol
    md1tau = md1.rqsol[1,]
    md1tau_beta = rbind(md1tau, md1.rqsol[-(1:3),])
    md1tau_beta = t(md1tau_beta)
  } else{
    # the conditional quantile function is nonlinear
    list_nlrq_1_model <- vec_nlrq_1_tau <- NULL

    # take the parallel computing if possible
    if("cluster" %in% class(clnodes)) {
      list_1 =  parLapply(cl=clnodes, X=tauvec,
                          fun=function(tau, moCondQuant_1, data1,start_1){
                            md1.nls.i = nlrq(formula = moCondQuant_1,
                                             data=data1,
                                             start = start_1, tau = tau)},
                          moCondQuant_1,
                          data1, start_1)
      for(i in 1:length(tauvec)){
        list_nlrq_1_model[[length(list_nlrq_1_model)+1]]<- list_1[[i]]
        vec_nlrq_1_tau <- c(vec_nlrq_1_tau, tauvec[i])
      }
    } else {
      # no parallel
      for(i in 1:length(tauvec)){
        md1.nls.i = nlrq(formula = moCondQuant_1, data=data1,
                         start = start_1, tau = tauvec[i])
        list_nlrq_1_model[[length(list_nlrq_1_model)+1]]<- md1.nls.i
        vec_nlrq_1_tau <- c(vec_nlrq_1_tau, tauvec[i])
      }
    }
  }

  # 2. to generate pseudo response ####
  # need random quantiles in [0,1]
  nobs <- nrow(raw.data)
  u0 <- runif(nobs, min = 0.002, max=0.998)
  u1 <- runif(nobs, min = 0.002, max=0.998)

  # According to whether the model is nonlinear, select the parameters
  # check beta_nlrq or beta_rq is applicable here
  # y.a.0 and y.a.1 are of the same length as the number off observations
  if(!nlCondQuant_0){
    betahat_index_0 <- unlist(lapply(u0, Betahat_rq, sol = md0tau_beta))
    mo_0_str <- str_replace_all(Reduce(paste, deparse(moCondQuant_0)), " ","")
    mo_0_terms <- strsplit(mo_0_str, "[~]")[[1]][2]
    responseNm <- strsplit(mo_0_str, "[~]")[[1]][1]
    dt0 <- model.frame(moCondQuant_0, raw.data)
    dt0 <- cbind(Intercept=1,dt0[,!(names(dt0) %in% c(responseNm))])
    beta0 <- md0tau_beta[betahat_index_0, 2:ncol(md0tau_beta), drop=FALSE]
    y.a.0 <- apply(dt0 * beta0, 1, sum)
    names(y.a.0)<-NULL
  } else {
    betahat_index_0 <- unlist(lapply(u0, Betahat_nlrq, Taus = vec_nlrq_0_tau))
    y.a.0 <- sapply(seq(nrow(raw.data)), FUN = function(i) {
      predict(list_nlrq_0_model[[ betahat_index_0[i]]], newdata = raw.data[i,])
    })
  }

  if(!nlCondQuant_1){
    betahat_index_1 <- unlist(lapply(u1, Betahat_rq, sol = md1tau_beta))
    mo_1_str <- str_replace_all(Reduce(paste, deparse(moCondQuant_1)), " ","")
    mo_1_terms <- strsplit(mo_1_str, "[~]")[[1]][2]
    responseNm <- strsplit(mo_1_str, "[~]")[[1]][1]
    dt1 <- model.frame(moCondQuant_1, raw.data)
    dt1 <- cbind(Intercept=1,dt1[,!(names(dt1) %in% c(responseNm))])
    beta1 <- md1tau_beta[betahat_index_1, 2:ncol(md1tau_beta), drop=FALSE]
    y.a.1 <- apply(dt1 * beta1, 1, sum)
    names(y.a.1)<-NULL
  } else {
    betahat_index_1 <- unlist(lapply(u1, Betahat_nlrq, Taus = vec_nlrq_1_tau))
    y.a.1 <- sapply(seq(nrow(raw.data)), FUN = function(i) {
      predict(list_nlrq_1_model[[betahat_index_1[i]]], newdata = raw.data[i,])
    })
  }

  return (list(y.a.0 = y.a.0,
               y.a.1 = y.a.1,
               nlCondQuant_0=nlCondQuant_0,
               nlCondQuant_1=nlCondQuant_1))
}



Betahat_rq<- function(tau, sol){ # The solution b(tau_i)
  #prevails from tau_i to tau_i+1
  Taus <-sol[,1]
  inds <- Taus<tau
  which.max(Taus[inds])
}

Betahat_nlrq<- function(tau, Taus){
  # The solution b(tau_i)
  #prevails from tau_i to tau_i+1
  ind <-  which.min(abs(Taus-tau))
  return(ind)
}



Test <- function(arg1,tau,qrwts.y,qrwts.a,y,y.a){
  Check_func<- function(tau,u){ u *(tau-as.numeric(u<0))  }

  return(sum(qrwts.y * Check_func(tau,(y - arg1)))
         +sum(qrwts.a * Check_func(tau,(y.a - arg1)) ))
}


#' @title  The Doubly Robust Quantile Estimator for a Given Treatment Regime
#' @description Given a fixed treatment regime, this doubly robust estimator
#' estimates the marginal quantile of responses when it is followed by
#' every unit in the target population. It took advantages of conditional
#' quantile functions for different treatment levels when they are available.
#' 
#' @inheritParams abso_diff_est
#' 
#' @param tau The quantile of interest
#' @param y.a.0 Estimated conditional potential outcome given that treatment = 0,
#' which can be calculated by the function \code{augX}.
#' @param y.a.1 Estimated conditional potential outcome given that treatment = 1,
#' which can be calculated by the function \code{augX}.
#' @param num_min logical. If \code{TRUE}, the number of global minimizers for the 
#'                objective function is returned.
#' @export
#' @details
#' The double robustness property means that it can consistently estimate
#' the marginal quantile when either the propensity score model is correctly
#' specified, or the conditional quantile function is correctly specified.
#' @seealso \code{\link{augX}}

dr_quant_est <- function(beta, x, y, a, prob, tau, y.a.0, y.a.1, num_min =FALSE){
  g <- as.numeric( x %*% beta > 0)
  c_consis <- a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  qrwts.y <- c_consis*wts # wts for the Y term
  qrwts.a <- 1- qrwts.y # wts for the Y^{\ast\ast} term
  y.a = g*y.a.1 + (1-g)*y.a.0
  candidate = c(y,y.a)

  # Qtilde is a vector, indicating value of the function Qtilde at every candidate point.
  # Since function Qtilde is piecewise linear,
  # global minimum can be attainde in the set "candidate".
  Qtilde = sapply(candidate, Test, tau, qrwts.y, qrwts.a, y, y.a)
  if(num_min)
    return(sum(Qtilde==min(Qtilde)))
  else
    return (candidate[which.min(Qtilde)])
}
