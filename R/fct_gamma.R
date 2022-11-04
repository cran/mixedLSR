#' Internal Posterior Calculation
#'
#' @inheritParams fct_em
#' @param pi_vec A vector of mixing probabilities for each cluster label.
#' @param N The sample size.
#' @param alpha A positive constant DPP parameter.
#' @param beta A positive constant DPP parameter.
#' @param y_sparse Should Y coefficients be treated as sparse?
#' @param rank The rank, if known.
#' @param max_rank The maximum allowed rank.
#'
#' @keywords internal
#'
#' @importFrom purrr safely
#'
#' @return A list with the posterior, coefficients, and estimated covariance.
#'
fct_gamma <- function(x, y, k, N, clust_assign, pi_vec, lambda, alpha, beta, y_sparse, rank, max_rank){

  s_mean <- purrr::safely(mean)
  safe_dpp <- purrr::safely(fct_dpp)
  safe_rank <- purrr::safely(fct_rank)
  p <- dim(x)[2]
  m <- dim(y)[2]
  val_frac <- 0.2
  grid <- seq(-5,4,1)
  if(length(lambda)==1){
    lambda <- rep(lambda, k)
  }

  gamma <- NULL
  A <- NULL
  sig_vec <- NULL
  for (i in seq(1,k)){

    cluster_rows <- which((clust_assign==i))
    n_k <- length(cluster_rows)
    x_k <- x[cluster_rows,]
    y_k <- y[cluster_rows, ]
    lambda_k <- lambda[i]
    eta_k <- 3


    if (n_k > 1){


      sigma_hat <- fct_sigma(y_k, n_k, m)

      if(is.null(rank)){
        rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta_k)
        rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
        rank_hat <- min(rank_hat, max_rank)
      } else {
        rank_hat <- rank
      }

      model_k_att <- safe_dpp(y_k, x_k, rank_hat, lambda_k, alpha, beta, sigma_hat, "grLasso", y_sparse)
      if(is.null(model_k_att$error)){
        model_k <- model_k_att$result

        A_k <- model_k$Ahat
        sigvec <- model_k$sigvec
        mu_mat <- cbind(x,1) %*% A_k
        gam <- fct_log_lik(mu_mat, sigvec, y, N, m) + log(pi_vec[i])

        gamma <- cbind(gamma, gam)
        A <- c(A,list(A_k))
        sig_vec <- c(sig_vec,list(sigvec))
      } else {
        gamma <- cbind(gamma, rep(-Inf,N))
        A <- c(A,list(NULL))
        sig_vec <- c(sig_vec,list(NULL))
      }

    } else {

      gamma <- cbind(gamma, rep(-Inf,N))
      A <- c(A,list(NULL))
      sig_vec <- c(sig_vec,list(NULL))

    }

  }
  return(list(gamma = gamma, A = A, sig_vec = sig_vec))
}
