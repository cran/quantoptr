#' Quantile- and Mean-Optimal Treatment Regimes
#' 
#' @description 
#'  \code{quantoptr} implements a variety of algorithms for estimating
#'  one-stage optimal treatment regimes(TRs)
#'  from either randomized controlled studies or observational studies
#'   under the quantile, the mean, and the mean absolute difference(MAD) criteria, respectively.
#'   It also considers the problem of estimating dynamic quantile- and mean-optimal TRs when 
#'   the decision-making process has two stages.
#'  
#'  
#'   
#' @details    
#' Treatment regimes (TRs) are decision rules that recommend treatments based on 
#' individual's characteristics and medical history. The problem of learning a policy
#' from either randomized controlled studies or observational studies is an active
#' research area with potential applications in both social and medical science. 
#' 
#' Depending on the objective, the shape of the probability distribution of outcome,  
#' and other restrictions of a particular application, a practitioner who wish
#' to construct an optimal individualized treatment regime can choose
#'  from a variety of optimality criteria. We assume that  larger values of the
#'   outcome variable are more favorable. Mean-optimal TR maximizes the average
#'   outcome in the potential population; quantile-optimal TR maximizes the marginal
#'   quantile in the potential population; and mean absolute difference-optimal TR
#'   minimizes the MAD, a measurement of statistical dispersion.
#'  
#'  The \code{quantoptr} package focuses on a class of direct-search estimators for
#'  the aforementioned criteria.
#'  Unlike regression-based methods, such as Q-learning (Murphy 2005), there is no need for an outcome model. 
#'  Rather, the direct-search estimators cast the problem as a missing data problem 
#'  and applies optimization directly on a prespecified class of rules. 
#'  More specifically, for one-stage problems, this package provides estimators for quantile-, mean-, 
#'  and MAD-optimal treatment regimes (Wang et. al. 2016, Zhang et.al. 2012, 2013). 
#'  Also, it provides a doubly robust estimator for estimating the
#'  quantile-optimal TR based on conditional quantile regression functions (Wang et. al. 2016).
#'  For two-stage problems, this package provides estimators for quantile- and 
#'  mean-optimal treatment regimes (not the doubly robust version). 
#'
#'  The functions that directly estimate optimal TRs all begin with capital letters,
#'  while other supporting functions begin with lower-case letters.
#' 
#' 
#' 
#' @references 
#' \insertRef{zhang2012robust}{quantoptr}
#' 
#' \insertRef{zhang2013robust}{quantoptr}
#' 
#' \insertRef{wang2016quant}{quantoptr}
#' 
#' \insertRef{murphy2005generalization}{quantoptr}
#' 
"_PACKAGE"
#> [1] "_PACKAGE"
