#' Internal Double Penalized Projection Function
#'
#' @inheritParams fct_gamma
#' @param sigma An estimated standard deviation
#' @param ptype A group penalized regression penalty type. See \link[grpreg]{grpreg}.
#'
#' @keywords internal
#'
#' @importFrom grpreg grpreg cv.grpreg
#' @importFrom stats sd
#'
#' @return A list containing estimated coefficients, covariance, and penalty parameters.
#'
fct_dpp <- function(y, x, rank, lambda = NULL, alpha = 2*sqrt(3), beta = 1, sigma, ptype = "grLasso", y_sparse = TRUE){
  n <- dim(x)[1]
  p <- dim(x)[2]
  m <- dim(y)[2]
  group  <- rep(seq(1,(p+1)),rank)
  Y_thresh <- matrix(0, ncol = m, nrow = n)
  lambda_store <- c(0,0)


  thresh_1 <- sigma^2*(n+alpha*sqrt(n*log(max(p,m))))
  j0 <- which(apply(y, 2, function(x){sum(x^2)}) >= thresh_1)
  Y0 <- Y_thresh
  if(y_sparse){
    Y0[,j0] <- y[,j0]
  } else {
    Y0 <- y
  }

  V0  <- svd(Y0,nu=rank,nv=rank)$v

  XX  <- kronecker(diag(rep(1,rank)),cbind(x,1))
  YY  <- y %*% V0
  YY  <- as.vector(YY)
  if(is.null(lambda )){
    fit1_cv <- grpreg::cv.grpreg(XX,YY,group, penalty= ptype, nfolds = 10, family="gaussian")
    fit1  <- grpreg::grpreg(XX,YY,group,lambda=fit1_cv$lambda.min, penalty= ptype, family="gaussian")
    lambda_store[1] <- fit1_cv$lambda.min
  } else {
    fit1  <- grpreg::grpreg(XX,YY,group,lambda=lambda , penalty= ptype, family="gaussian")
    lambda_store[1] <- lambda
  }

  B1 <- matrix(fit1$beta[-1],nrow=p+1,ncol=rank)
  B1[p+1,] <-B1[p+1,]+fit1$beta[1]
  XB <-cbind(x,1) %*% matrix(B1,nrow=p+1,ncol=rank)
  U1 <- svd(XB,nu=rank,nv=rank)$u

  thresh_2 <- beta*sigma^2*(rank + 2*sqrt(3*rank*log(max(p,m))) + 6*log(max(p,m)))
  j1_tmp <- which(apply(y, 2, function(x){sum((t(U1)%*%matrix(x))^2)}) > thresh_2)
  j1 <- sort(unique(c(j0, j1_tmp)))
  Y1 <- Y_thresh

  if(y_sparse){
    Y1[,j1] <- y[,j1]
  } else {
    Y1 <- y
  }


  tmp <- U1%*%t(U1)%*%Y1
  V1 <- svd(tmp,nu=rank,nv=rank)$v

  YY  <- y %*% V1
  YY  <- as.vector(YY)

  if(is.null(lambda )){
    fit2_cv <- grpreg::cv.grpreg(XX,YY,group, penalty= ptype, family="gaussian")
    fit2  <- grpreg::grpreg(XX,YY,group,lambda=fit2_cv$lambda.min, penalty= ptype, nfolds = 10, family="gaussian")
    lambda_store[2] <- fit2_cv$lambda.min
  } else {
    fit2  <- grpreg::grpreg(XX,YY,group,lambda=lambda , penalty= ptype, nfolds = 10, family="gaussian")
    lambda_store[2] <- lambda
  }
  B2 <- matrix(fit2$beta[-1],nrow=p+1,ncol=rank)
  B2[p+1,] <-B2[p+1,]+fit2$beta[1]

  Ahat <- B2 %*% t(V1)
  sigvec  <- apply(y- cbind(x,rep(1,nrow(x)))%*% Ahat, 2, stats::sd)
  return(list(Ahat=Ahat,sigvec=sigvec,lambda_store=lambda_store))
}
