#' Internal EM Algorithm
#'
#' @inheritParams fct_alt_optimize
#' @param lik_track A vector storing the log-likelihood by iteration.
#'
#' @keywords internal
#'
#' @return A mixedLSR model.
#'
fct_em <- function(x, y, k, lambda, clust_assign, lik_track, em_iter, verbose){
  iter <- 0
  N <- nrow(x)
  changed <- Inf
  while(changed > 0 & iter < em_iter){
    iter <- iter + 1

    pi_vec <- fct_pi_vec(clust_assign, k, N)
    gamma_model <- fct_gamma(x = x, y = y, k = k, N = N, clust_assign = clust_assign,
                             pi_vec = pi_vec, lambda = lambda, alpha = 2*sqrt(3),
                             beta = 1, y_sparse = TRUE, rank = NULL, max_rank = 3)
    gamma <- gamma_model$gamma
    A <- gamma_model$A
    sig_vec <- gamma_model$sig_vec
    weighted_ll <- fct_weighted_ll(gamma)
    lik_track <- rbind(lik_track,data.frame(iter = lik_track$iter[nrow(lik_track)]+1,ll=weighted_ll,type="EM"))


    clust_assign_old <- clust_assign
    clust_assign <- apply(gamma,1,which.max)
    old_ll <- weighted_ll


    changed <- sum(clust_assign != clust_assign_old)
    if(verbose){cat(".")}
    if(changed > 0){
      lambda_old <- lambda
      lambda <- fct_select_lambda(x, y, k, clust_assign, initial = FALSE, verbose = verbose)
      empty_lam <- which(lambda==0)
      lambda[empty_lam] <- lambda_old[empty_lam]

    }

  }
  if(verbose){cat("\n")}
  return(list(assign = clust_assign, lambda = lambda, lik_track = lik_track,
              weighted_ll = weighted_ll))

}
