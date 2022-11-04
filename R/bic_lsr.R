#' Compute Bayesian information criterion for a mixedLSR model
#'
#' @param a A list of coefficient matrices.
#' @param n The sample size.
#' @param llik The log-likelihood of the model.
#'
#'
#' @return The BIC.
#' @export
#'
#' @examples
#' n <- 50
#' simulate <- simulate_lsr(n)
#' model <- mixed_lsr(simulate$x, simulate$y, k = 2, init_lambda = c(1,1), alt_iter = 0)
#' bic_lsr(model$A, n = n, model$llik)
bic_lsr <- function(a, n, llik){
  k <- length(a)
  nu <- 0
  p <- nrow(a[[1]])*ncol(a[[1]])
  for(i in seq(1,k)){
    j <- sum(a[[i]]!=0)
    nu <- nu + j
  }
  -2*llik+nu*log(n)
}


