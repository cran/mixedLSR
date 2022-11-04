#' Internal Simulated Annealing Function
#'
#' @inheritParams fct_alt_optimize
#' @param init_assign An initial clustering assignment.
#' @param track A likelihood tracking vector.
#'
#' @keywords internal
#'
#' @return An updated clustering vector.
#'
fct_sim_anneal <- function(x, y, k, init_assign, lambda, temp, mu, eps, accept_prob, sim_N, track, anneal_iter = 1e3, verbose){
  total_iter <- 0
  count <- 0

  a_t <- a_b <- a_c <- init_assign
  j_b <- j_c <- fct_j_lik(x, y, k, init_assign, lambda)
  t <- temp


  while(t >= eps & total_iter <= anneal_iter){
    total_iter <- total_iter + 1
    a_t <- fct_new_assign(a_b, k, accept_prob)
    j_t <- fct_j_lik(x, y, k, a_t, lambda)

    if(j_t <= j_c){
      a_c <- a_t
      j_c <- j_t

      if(j_t >= j_b){
        count <- count + 1
      } else {
        changed <- sum(a_b != a_t)
        j_b <- j_t
        a_b <- a_t
        count <- 0
        total_iter <- 0

      }
    } else {
      u <- runif(1)
      b <- exp(-(j_t - j_c)/t)
      if(u <= b){
        j_c <- j_t
        a_c <- a_t
      }
    }
    if(count >= sim_N){
      t <- mu*t
    }
    track <- rbind(track,data.frame(iter=track$iter[nrow(track)]+1, ll=-j_b, type = "sim"))
    if((total_iter %% 50 == 0) & verbose){cat(".")}
  }
  if(verbose){cat("\n")}

  if(sum(a_b != init_assign) > 0){

    lambda_old <- lambda
    lambda <- fct_select_lambda(x, y, k, a_b, initial = FALSE, verbose = verbose)
    empty_lam <- which(lambda==0)
    lambda[empty_lam] <- lambda_old[empty_lam]
  }
  return(list(assign = a_b, lik_track = track, weighted_ll = -j_b, lambda = lambda))
}
