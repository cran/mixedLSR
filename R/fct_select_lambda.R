#' Internal Penalty Parameter Selection Function.
#'
#' @inheritParams fct_alt_optimize
#' @param initial An initial penalty parameter.
#' @param type A type.
#'
#' @keywords internal
#'
#' @importFrom purrr safely
#' @importFrom stats median
#'
#' @return A selected penalty parameter.
#'
fct_select_lambda <- function(x, y, k, clust_assign = NULL, initial = FALSE, type = "all", verbose){
  max_rank <- 3
  safe_rank <- purrr::safely(fct_rank)
  if(verbose & initial){cat("Selecting Lambda")}
  if(initial){
    M <- 50
    clust_assign <- fct_initialize(k,nrow(x))
  } else{
    M <- 2
  }
  store <- array(0, dim = c(M,k,2))

  for(i in seq(1,M)){
    if(initial){
      clust_assign <- fct_initialize(k,nrow(x))
    }
    if(verbose & initial){cat(".")}
    for(j in seq(1,k)){
      if(length(which(clust_assign==j))>2){
        x_k <- x[which(clust_assign==j),]
        y_k <- y[which(clust_assign==j),]
        sigma_hat <- fct_sigma(y_k, nrow(y_k), ncol(y_k))

        rank_sest <- safe_rank(x_k, y_k, sigma_hat, eta = 3)
        rank_hat <- ifelse(is.null(rank_sest$result),1,rank_sest$result)
        rank_hat <- min(rank_hat, max_rank)

        store[i,j,] <- fct_dpp(y_k, x_k, rank = rank_hat, lambda = NULL,
                                 alpha = 2*sqrt(3), beta = 1, sigma = sigma_hat,
                                 ptype = "grLasso", y_sparse = TRUE)$lambda_store
      }

    }
  }
  if(verbose & initial){cat("\n")}
  store_mat <- rbind(store[,,1],store[,,2])
  if(type == "single"){
    lambda <- stats::median(store_mat)
  } else {
    lambda <- apply(store_mat, 2, stats::median)
  }
  return(lambda)
}
