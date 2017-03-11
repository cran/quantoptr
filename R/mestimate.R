#' @title The Mean-Optimal Treatment Regime Wrapper Function
#' @description
#' The wrapper function for mean-optimal treatment regime that calls a genetic algorithm.
#' This function supports the \code{\link{IPWE_Mopt}} function.
#' 
#' @param x    a matrix of observed covariates from the sample. 
#' Notice that we assumed the class of treatment regimes is linear.
#' @param nvars an integer. The number of parameters indexing a treatment regime.
#' @param hard_limit  logical. This logical variable determines if the
#' max.generations variable is a binding constraint for genoud.
#' @param max logical. If \code{max=TRUE}, it indicates we wish to maximize the marginal
#' mean; If \code{max=FALSE}, we wish to minimize the marginal mean. The default is \code{TRUE}.
#' 
#' @inheritParams abso_diff_est
#' @inheritParams IPWE_Qopt
#'
#'
#' @seealso The function \code{\link{IPWE_Mopt}} is based on this function.
#' @importFrom  rgenoud genoud
#' @import stats
#' 
#' @references 
#' \insertRef{zhang2012robust}{quantoptr}
#' 
#' @export
mestimate<-function(x,y,a,prob,p_level,nvars,
                    hard_limit=FALSE,max=TRUE,
                    cl.setup = 1, s.tol=0.0001, it.num = 8, pop.size = 3000)
{
  # cl.setup is the number of cores to use
  if(cl.setup>1){
    # parallel computing option
    if(!(get_os()== "windows"))
      clnodes <- parallel::makeCluster(cl.setup, type="PSOCK")
    else if((get_os()== "osx") | (get_os()== "linux"))
      clnodes <- parallel::makeForkCluster(nnodes = getOption("mc.cores", cl.setup))
    else
      # no parallel
      clnodes <- FALSE
  } else {
    # no parallel
    clnodes <- FALSE
  }

  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est <- genoud(fn=mean_est,nvars=nvars,x=x,y=y,a=a,prob=prob,
                print.level=p_level,max=max,
                pop.size=pop.size,
                wait.generations=it.num,
                gradient.check=FALSE, BFGS=FALSE,
                P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
                Domains=Domains,
                starting.values=rep(0,nvars),
                hard.generation.limit=hard_limit,
                solution.tolerance=s.tol,
                optim.method="Nelder-Mead",
                cluster = clnodes)

  if("cluster" %in%class(clnodes)) { parallel::stopCluster(clnodes) }

  #########  estimated  coefficient ####################
  coefficient <- est$par
  coefficient <- scalar1(coefficient)

  output.v<-list(coefficient=coefficient,hatM = est$value)
  return(output.v)
}
