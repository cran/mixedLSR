#' Simulate Heterogeneous, Low-Rank, and Sparse Data
#'
#' @param N The sample size, default = 100.
#' @param k The number of groups, default = 2.
#' @param p The number of predictor features, default = 30.
#' @param m The number of response features, default = 35.
#' @param b The signal-to-noise ratio, default = 1.
#' @param d The singular value, default = 20.
#' @param h The lower bound for the singular matrix simulation, default = 0.2.
#' @param case The covariance case, "independent" or "dependent", default = "independent".
#'
#' @importFrom stats rnorm
#'
#' @return A list of simulation values, including x matrix, y matrix, coefficients and true clustering assignments.
#' @export
#'
#' @examples
#' simulate_lsr()
simulate_lsr <- function(N = 100, k = 2, p = 30, m = 35, b = 1, d = 20, h = 0.2, case = "independent"){
  rank <- 1
  prob <- rep(1/k,k)
  int <- cumsum(prob)
  rand_assign <- stats::runif(N)

  clust_assign_true <- sapply(rand_assign, function(x){min(which(x <= int))})
  n <- sapply(seq(1,k), function(x){sum(clust_assign_true==x)})

  x <- matrix(0, nrow = N, ncol = p)
  y <- matrix(0, nrow = N, ncol = m)
  a <- NULL
  comp <- rep(1:3, length.out = k)
  for(i in seq(1, k)){

    u1 <- c(sample(c(1,-1), 5, replace = TRUE)*runif(5,h,1),rep(0,20))
    u1 <- u1/norm(u1, type = "2")
    u2 <- c(rep(0,5),sample(c(1,-1), 5, replace = TRUE)*runif(5,h,1),rep(0,15))
    u2 <- u2/norm(u2, type = "2")
    u3 <- c(rep(0,10),sample(c(1,-1), 5, replace = TRUE)*runif(5,h,1),rep(0,10))
    u3 <- u3/norm(u3, type = "2")
    u <- matrix(cbind(u1,u2,u3)[,comp[i]],ncol = rank)

    v1 <- c(sample(c(1,-1), 10, replace = TRUE),rep(0,15))
    v1 <- v1/norm(v1, type = "2")
    v2 <- c(rep(0,12),sample(c(1,-1), 10, replace = TRUE),rep(0,3))
    v2 <- v2/norm(v2, type = "2")
    v3 <- c(sample(c(1,-1), 6, replace = TRUE),v1[7:8],-v1[9:10],sample(c(1,-1), 2, replace = TRUE),
            -v1[13:14],v1[15:16],rep(0,9))
    v3 <- v3/norm(v3, type = "2")
    v <- matrix(cbind(v1,v2,v3)[,comp[i]],ncol = rank)

    if(p > 25){u <- rbind(u,matrix(0,nrow = p-25,ncol=rank))}
    if(m > 25){v <- rbind(v,matrix(0,nrow = m-25,ncol=rank))}

    u_shuffle <- sample(1:p, p)
    v_shuffle <- sample(1:m, m)

    u <- matrix(u[u_shuffle,], ncol = rank)
    v <- matrix(v[v_shuffle,], ncol = rank)

    C <- u%*%d%*%t(v)

    if(case=="independent"){rho = 0} else {rho = sample(seq(0.125,0.625,0.125),1)}

    Sigma <- matrix(1,p,p)
    for(j in 1:p) for (k in 1:p) Sigma[j,k] <- rho^abs(j-k)


    X <- MASS::mvrnorm(n[i],rep(0,p),Sigma)
    E <- matrix(stats::rnorm(n[i]*m, sd = sum(diag(t(C)%*%Sigma%*%C))/(n[i]*m*b)),n[i],m)
    Y <- X%*%C + E

    a <- c(a,list(C))
    x[clust_assign_true==i,] <- X
    y[clust_assign_true==i,] <- Y
  }

  list(x = x, y = y, a = a, true = clust_assign_true)
}
