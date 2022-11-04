#' Internal Alternating Optimization Function
#'
#' @inheritParams mixed_lsr
#' @param clust_assign The current clustering assignment.
#' @param lambda A vector of penalization parameters.
#'
#' @keywords internal
#'
#' @return A final fit of mixedLSR
#'
fct_alt_optimize <- function(x, y, k, clust_assign, lambda, alt_iter, anneal_iter, em_iter, temp, mu, eps, accept_prob, sim_N, verbose){
  N <- nrow(x)
  p <- ncol(x)
  m <- ncol(y)

  if(is.null(clust_assign)){
    clust_assign <- fct_initialize(k,N)
  }
  iter <- 0
  lik_track <- data.frame(iter=0, ll = -Inf, type = "EM")
  changed <- Inf

  # Select Initial Lambda
  if(is.null(lambda)){
    lambda <- fct_select_lambda(x, y, k, clust_assign = NULL, initial = TRUE, verbose = verbose)
  }

  while(iter < alt_iter){

    # EM
    if(verbose){cat("EM Step")}
    clust_assign_old <- clust_assign
    model_em <- fct_em(x = x, y = y, k = k, lambda = lambda, clust_assign = clust_assign,
                       lik_track = lik_track, em_iter = em_iter, verbose)
    clust_assign <- model_em$assign
    lambda <- model_em$lambda
    lik_track <- model_em$lik_track
    weighted_ll <- model_em$weighted_ll


    # Simulated Annealing
    if(verbose){cat("Simulated Annealing Step")}
    model_anneal <- fct_sim_anneal(x = x, y = y, k, init_assign = clust_assign,
                                   lambda = lambda, temp = temp, mu = mu, eps = eps,
                                   accept_prob = accept_prob, sim_N = sim_N, track = lik_track,
                                   anneal_iter = anneal_iter, verbose)
    clust_assign <- model_anneal$assign
    lambda <- model_anneal$lambda
    lik_track <- model_anneal$lik_track
    weighted_ll <- model_anneal$weighted_ll

    changed <- sum(clust_assign != clust_assign_old)
    iter <- iter + 1
    if(verbose){cat(paste("Full Cycle", iter,"\n"))}


  }

  # final fit
  if(verbose){cat("Computing Final Model...\n")}
  pi_vec <- fct_pi_vec(clust_assign, k, N)
  final_model <- fct_gamma(x, y, k, N, clust_assign, pi_vec, lambda,
                           alpha = 2*sqrt(3), beta = 1, y_sparse = TRUE,
                           max_rank = 3, rank = NULL)
  gamma <- final_model$gamma
  A <- final_model$A
  sig_vec <- final_model$sig_vec
  weighted_ll <- fct_weighted_ll(gamma)

  return(list(ll = weighted_ll, assign = clust_assign, A = A, sig_vec = sig_vec,
              ll_store = lik_track, iter = iter))

}
