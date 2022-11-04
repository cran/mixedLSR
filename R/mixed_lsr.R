#' Mixed Low-Rank and Sparse Multivariate Regression for High-Dimensional Data
#'
#' @param x A matrix of predictors.
#' @param y A matrix of responses.
#' @param k The number of groups.
#' @param nstart The number of random initializations, the result with the maximum likelihood is returned.
#' @param init_assign A vector of initial assignments, NULL by default.
#' @param init_lambda A vector with the values to initialize the penalization parameter for each group, e.g., c(1,1,1). Set to NULL by default.
#' @param alt_iter The maximum number of times to alternate between the classification expectation maximization algorithm and the simulated annealing algorithm.
#' @param anneal_iter The maximum number of simulated annealing iterations.
#' @param em_iter The maximum number of EM iterations.
#' @param temp The initial simulated annealing temperature, temp > 0.
#' @param mu The simulated annealing decrease temperature fraction. Once the best configuration cannot be improved, reduce the temperature to (mu)T, 0 < mu < 1.
#' @param eps The final simulated annealing temperature, eps > 0.
#' @param accept_prob The simulated annealing probability of accepting a new assignment 0 < accept_prob < 1. When closer to 1, trial assignments will only be small perturbation of the current assignment. When closer to 0, trial assignments are closer to random.
#' @param sim_N The simulated annealing number of iterations for reaching equilibrium.
#' @param verbose A boolean indicating whether to print to screen.
#'
#' @return A list containing the likelihood, the partition, the coefficient matrices, and the BIC.
#' @export
#'
#' @examples
#' simulate <- simulate_lsr(50)
#' mixed_lsr(simulate$x, simulate$y, k = 2, init_lambda = c(1,1), alt_iter = 0)
mixed_lsr <- function(x, y, k, nstart = 1, init_assign = NULL, init_lambda = NULL, alt_iter = 5, anneal_iter = 1e3,
                   em_iter = 1e3, temp = 1e3, mu = 0.95, eps = 1e-6, accept_prob = 0.95, sim_N = 200, verbose = TRUE){

  N <- nrow(x)
  assignments <- NULL
  likelihood <- NULL
  ll_store <- NULL
  A <- NULL
  for (i in seq(1,nstart)){
    if(verbose){cat(paste("mixedLSR Start:",i,"\n"))}

    model <- fct_alt_optimize(x, y, k, init_assign, init_lambda, alt_iter, anneal_iter, em_iter, temp, mu, eps, accept_prob, sim_N, verbose)
    likelihood <- c(likelihood,model$ll)
    assignments <- c(assignments,list(model$assign))
    ll_store <- c(ll_store,list(model$ll_store))
    A <- c(A,list(model$A))

  }

  best <- which.max(likelihood)
  best_bic <- bic_lsr(A[[best]], N, likelihood[best])
  result <- list(llik = likelihood[best], assign = assignments[[best]], a = A[[best]], BIC = best_bic)
  if(verbose){cat("Done! \n")}
  return(result)
}
