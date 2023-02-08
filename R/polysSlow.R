#' @importFrom minqa bobyqa
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom stats optimize
#' @importFrom stats cor
polysSlow <- function(x, M, w, ML=FALSE) {
  # polyserial log likelihood function
  polysLnL <- function(x,M,rho,theta,w) {
    R <- (1-rho^2)^0.5
    Qp2 <- (theta[M+2] - rho*x) / R
    Qp1 <- (theta[M+1] - rho*x) / R
    sum(w * dnorm(x, log=TRUE)) + sum(w * log(pnorm(Qp2) - pnorm(Qp1)))
  }
  # weighted normalization of x
  fixx <- function(x,w) {
    mux <- sum(x*w)
    sdx <- sum(w*(x-mux)^2)
    (x-mux)/sqrt(sdx)
  }
  # to allow the correltion to be fit on R instead of [0,1] use this remapping
  mapCor <- function(v) {
    tanh(v)
  }
  # puts the theta into the correct format for the olysLnL function
  mapTheta <- function(v) {
    vv <- cumsum(c(v[1],exp(v[-1])))
    c(NA, -Inf, vv, Inf)
  }
  
  optF <- function(x,M,w) {
    w <- w/sum(w)
    fx <- fixx(x,w)
    function(par) {
      res <- polysLnL(fx, M, mapCor(par[1]), mapTheta(par[-1]), w)
      ifelse(res==-Inf, .Machine$double.xmax, -1*res)
    }
  }

  optFc <- function(x,M,w,theta0) {
    w <- w/sum(w)
    fx <- fixx(x,w)
    ftheta0 <- mapTheta(theta0)
    function(par) {
      res <- polysLnL(fx, M, mapCor(par[1]), ftheta0, w)
      ifelse(res==-Inf, .Machine$double.xmax, -1*res)
    }
  }
  
  imapTheta <- function(theta0) {
    c(theta0[1], log(theta0[-1] - theta0[-length(theta0)]))
  }
  # 
  imapCor <- function(cor) {
    atanh(cor)
  }
  # map cor so that the optimization can be over R
  mapCor <- function(v) {
    tanh(v)
  }

  M <- as.numeric(as.factor(M)) # make discrete values that are adjacent.
  uM <- sort(unique(M))
  theta0 <- sapply(uM[-length(uM)], function(z) qnorm( weighted.mean(M<=z, w) ) )
  if(ML) {
    bob <- bobyqa(par=c(imapCor(cor(x,M)), imapTheta(theta0)), fn=optF(x,M,w))
    return(mapCor(bob$par[1]))
  } else {
    opt <- optimize(optFc(x, M, w, imapTheta(theta0)), imapCor(cor(x,M)) + c(-3,3))
    return(mapCor(opt$minimum))
  }
}
