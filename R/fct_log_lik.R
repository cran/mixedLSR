#' Internal Log-Likelihood Function
#'
#' @param mu_mat The mean matrix.
#' @param sig_vec A vector of sigma.
#' @param y The output matrix.
#' @param N The sample size.
#' @param m The number of y features.
#'
#' @keywords internal
#'
#' @importFrom stats dnorm
#'
#' @return A posterior matrix.
#'
fct_log_lik <- function(mu_mat, sig_vec, y, N, m){

  gam <- array(dim = N)
  for(i in seq(1,N)){
    mu_i <- mu_mat[i,]
    gam[i] <- sum(log(sapply(seq(1,m), function(a){stats::dnorm(y[i,a], mean = mu_i[a], sd = sig_vec[a])})))
  }
  return(gam)
}
