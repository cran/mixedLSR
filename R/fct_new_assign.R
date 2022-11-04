#' Internal Perturb Function
#'
#' @param assign The current clustering assignments.
#' @param k The number of groups.
#' @param p The acceptance probability.
#'
#' @keywords internal
#'
#' @return A perturbed assignment.
#'
fct_new_assign <- function(assign, k, p){
  n <- length(assign)
  new_assign <- assign
  i <- 1
  flag <- FALSE

  while(!(i == n & flag)){
    if(i > n){i <- 1}
    u <- runif(1)
    if(u > p){
      flag <- TRUE
      cur_assign <- assign[i]
      new_assign[i] <- sample((seq(1,k))[-which(cur_assign==seq(1,k))],1)
    }
    i <- i + 1
  }
  return(new_assign)
}
