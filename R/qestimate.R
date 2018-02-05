#' @title The Quantile-Optimal Treatment Regime Wrapper Function
#' @description
#' The wrapper function for quantile-optimal treatment regime that calls a genetic algorithm.
#' This function supports the \code{\link{IPWE_Qopt}} function.
#' 
#' @param x    a matrix of observed covariates from the sample. 
#' Notice that we assumed the class of treatment regimes is linear.
#' @param tau a numeric value between 0 and 1. The quantile level of interest.
#' @param nvars an integer. The number of parameters indexing a treatment regime.
#' @param hard_limit  logical. This logical variable determines if the max.generations
#' variable is a binding constraint for \code{rgenoud::genoud()}.
#' 
#' 
#' @inheritParams abso_diff_est
#' @inheritParams IPWE_Qopt
#' 
#' 
#' @seealso The function \code{\link{IPWE_Qopt}} is based on this function.
#' 
#' 
#' @references 
#' \insertRef{wang2017quantile}{quantoptr}
#' 
#' @importFrom rgenoud genoud
#' @importFrom Rdpack insert_ref
#' @import stats
#' @export
qestimate<-function(tau, x, y, a, prob, p_level,nvars,hard_limit,max=TRUE,
                    cl.setup = 1,s.tol=0.0001, it.num = 8, pop.size = 3000)
{
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

  nvars <- ncol(x)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est <- genoud(fn=quant_est,
              nvars=nvars,
              x=x,
              y=y,a=a, prob=prob, tau=tau,
              print.level=p_level, 
              max=max,
              pop.size=pop.size,
              wait.generations=it.num,
              gradient.check=FALSE, BFGS=FALSE,
              P1=50, P2=50, P3=50, P4=50, P5=50,
              P6=50, P7=50, P8=50, P9=0,
              Domains=Domains,
              starting.values=rep(0,nvars),
              hard.generation.limit=hard_limit,
              solution.tolerance=s.tol,
              optim.method="Nelder-Mead",
              cluster = clnodes)

  if("cluster" %in%class(clnodes)) { parallel::stopCluster(clnodes) }

  #########  estimated  coefficient ##############
  coefficient<-est$par
  if (prod(coefficient==rep(0,nvars))==1)
    coefficient <- rep(0, nvars)
  else
    coefficient <- scalar1(coefficient)
  hatQ<-est$value

  output.v <- list(coefficient=coefficient, hatQ=hatQ)
  return(output.v)
}


