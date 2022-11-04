#' Internal Partition Initialization Function
#'
#' @param k The number of groups.
#' @param N The sample size.
#'
#' @keywords internal
#'
#' @importFrom stats runif
#'
#' @return A vector of assignments.
#'
fct_initialize <- function(k, N){

  init_int <- cumsum(rep(1/k,k))
  init_rand_assign <- stats::runif(N)
  sapply(init_rand_assign, function(x){min(which(x <=init_int))})

}
