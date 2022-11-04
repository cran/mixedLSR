#' Internal Rank Estimation Function
#'
#' @inheritParams fct_gamma
#' @param sigma An estimated noise level.
#' @param eta A rank selection parameter.
#'
#' @keywords internal
#'
#' @importFrom MASS ginv
#'
#' @return The estimated rank.
#'
fct_rank <- function(x, y, sigma, eta){
  m <- t(x)%*%x
  m_inv <- MASS::ginv(m)
  p <- x%*%m_inv%*%t(x)
  sv <- svd(p%*%y)$d
  max(which(sv >= sigma*eta))
}
