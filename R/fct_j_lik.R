#' Internal Likelihood Function
#'
#' @inheritParams fct_sim_anneal
#' @inheritParams fct_gamma
#' @param clust_assign A vector of cluster labels.
#'
#' @keywords internal
#'
#' @return The weighted log-likelihood
#'
fct_j_lik <- function(x, y, k, clust_assign, lambda, alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL){
  N <- nrow(x)
  pi_vec <- fct_pi_vec(clust_assign, k, N)
  gamma_model <- fct_gamma(x, y, k, N, clust_assign, pi_vec = pi_vec, lambda, alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE, max_rank = 3, rank = NULL)
  gamma <- gamma_model$gamma
  return(-fct_weighted_ll(gamma))
}
