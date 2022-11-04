#' Internal Pi Function
#'
#' @inheritParams fct_gamma
#'
#' @keywords internal
#'
#' @return A mixing vector.
#'
fct_pi_vec <- function(clust_assign, k, N){
  sapply(seq(1,k),function(x){sum(clust_assign==x)/length(clust_assign)})
}
