#' Internal Weighted Log Likelihood Function
#'
#' @param gamma A posterior matrix
#'
#' @keywords internal
#'
#' @return A weighted log likelihood vector
#'
fct_weighted_ll <- function(gamma){

  llik <- apply(gamma, 1, max)
  llik[is.infinite(llik)] <- min(llik[is.finite(llik)])
  sum(llik)

}
