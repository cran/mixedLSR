#' Internal Sigma Estimation Function
#'
#' @inheritParams fct_gamma
#' @param m The number of outcome variables.
#'
#' @keywords internal
#'
#' @importFrom stats median
#'
#' @return The estimated sigma.
#'
fct_sigma <- function(y, N, m){
  sv <- svd(y)$d
  stats::median(sv[sv!=0])/sqrt(max(N,m))
}
