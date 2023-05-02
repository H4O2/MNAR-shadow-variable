#Date: 2021/03/04
#### Doubly robust estimation of nonignorable normal data

simuTT <- function(N){
  X <- rnorm(N, mean=0, sd=1)
  Xs <- cbind(1,X) 
  Xr <- cbind(1,X)
  Xz <- cbind(1,X^2) 

  dr <- dim(Xr)[2]
  dz <- dim(Xz)[2]
  ds <- dim(Xs)[2]
  
  para = c(-0.5, 0.5, 0.4, 1, 0, 1, 0, -0.4, 1, 1)
  gamma <- para[1]
  alpha <- para[2:(ds+1)]
  beta1 <- para[(ds+2):(dr+ds+2)]
  beta2 <- para[(dr+ds+3):(dr+dz+ds+2)]
  beta12 <- beta1[1]
  sigmay2 <- para[dr+dz+ds+3]
  sigmaz2 <- para[dr+dz+ds+4]
  
  # Data-generating process with true propensity and true outcome model
  resp <- plogis(cbind(Xs,Xr,Xz)%*%c(alpha,gamma*beta1[-1],gamma*beta12*beta2) - gamma^2*(sigmay2+beta12^2*sigmaz2)/2)
  R <- rbinom(N,size=1,prob=resp)
  muZ <- Xz%*%beta2 + (R-1)*gamma*beta12*sigmaz2
  Z <- rnorm(N,mean=muZ,sd=sqrt(sigmaz2))
  muY <- cbind(Z,Xr)%*%beta1 + (R-1)*gamma*sigmay2
  Y <- rnorm(N, mean=muY, sd=sqrt(sigmay2))
  mu <- mean(Y)
  Y[R==0] <- NA
  
  data1 <- list(Y=Y,R=R,Z=Z,Xs=cbind(1,X),Xr=cbind(1,X),Xz=cbind(1,X^2))
  
  # IPW estimation 
  ipwpar <- optim(par = c(-0.03,-0.51,0.52,0.41),
                  fn = GMM,
                  f = IPWmrf, data=data1,
                  method = "BFGS", hessian = FALSE)$par
  # REG estimation
  regpar <- optim(par = c(-0.03,-0.51,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                 fn = GMM,
                 f = REGmrf, data=data1,
                 method = "BFGS", hessian = FALSE)$par
  # DR estimation
  drpar <- optim(par = c(-0.03,-0.51,0.52,0.41,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                 fn = GMM,
                 f = DRmrf, data=data1,
                 method = "BFGS", hessian = FALSE)$par
  # EFF estimation
  effres <- EFFmrf(drpar, data1)
  effpar <- apply(effres, 2, mean) + drpar[1:2]
  
  # MAR estimation
  lm1 <- lm(Y[R==1]~X[R==1])
  marmu <- mean(cbind(1,X)%*%lm1$coef)
  
  # Variance estimation
  ipwvar <- diag(VarEst(IPWmrf,ipwpar,data1))[1:2]
  regvar <- diag(VarEst(REGmrf,regpar,data1))[1:2]
  drvar <- diag(VarEst(DRmrf,drpar,data1))[1:2]
  effvar <- diag(VarEst(EFFmrf2,c(effpar,drpar),data1))[1:2]

  return(list(mu=mu, ipwmu=ipwpar[1], ipwor=ipwpar[2],
              regmu=regpar[1], regor=regpar[2],
              drmu=drpar[1], dror=drpar[2], 
              effmu=effpar[1], effor=effpar[2], 
              marmu=marmu,
              ipwvar=ipwvar, regvar=regvar, drvar=drvar, effvar=effvar)) 
}


simuFT <- function(N){
  X <- rnorm(N,mean=0,sd=1)
  Xr <- cbind(1,X)
  Xs <- cbind(1,X,X^2) 
  Xz <- cbind(1,X^2)
  
  dr <- dim(Xr)[2]
  dz <- dim(Xz)[2]
  ds <- dim(Xs)[2]
  
  para = c(-0.5, 0, 1, 1.5, 1, 0, 1, 0, -0.4, 1, 1)
  gamma <- para[1]
  alpha <- para[2:(ds+1)] 
  beta1 <- para[(ds+2):(dr+ds+2)]
  beta2 <- para[(dr+ds+3):(dr+dz+ds+2)]
  beta12 <- beta1[1]
  sigmay2 <- para[dr+dz+ds+3]
  sigmaz2 <- para[dr+dz+ds+4]
  
  # Data-generating process with false propensity and true outcome model
  resp <- plogis(cbind(Xs,Xr,Xz)%*%c(alpha,gamma*beta1[-1],gamma*beta12*beta2) - gamma^2*(sigmay2+beta12^2*sigmaz2)/2)
  R <- rbinom(N,size=1,prob=resp)
  muZ <- Xz%*%beta2 + (R-1)*gamma*beta12*sigmaz2
  Z <- rnorm(N,mean=muZ,sd=sqrt(sigmaz2))
  muY <- cbind(Z,Xr)%*%beta1 + (R-1)*gamma*sigmay2
  Y <- rnorm(N, mean=muY, sd=sqrt(sigmay2))
  mu <- mean(Y)
  Y[R==0] <- NA
  
  data1 <- list(Y=Y,R=R,Z=Z,Xs=cbind(1,X),Xr=cbind(1,X),Xz=cbind(1,X^2))
  
  # IPW estimation 
  ipwpar <- optim(par = c(-0.03,-0.51,0.52,0.41),
                  fn = GMM,
                  f = IPWmrf, data=data1,
                  method = "BFGS", hessian = FALSE)$par
  # REG estimation
  regpar <- optim(par = c(-0.03,-0.51,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                  fn = GMM,
                  f = REGmrf, data=data1,
                  method = "BFGS", hessian = FALSE)$par
  # DR estimation
  drpar <- optim(par = c(-0.03,-0.51,0.52,0.41,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                 fn = GMM,
                 f = DRmrf, data=data1,
                 method = "BFGS", hessian = FALSE)$par
  # EFF estimation
  effres <- EFFmrf(drpar, data1)
  effpar <- apply(effres, 2, mean) + drpar[1:2]
  
  # MAR estimation
  lm1 <- lm(Y[R==1]~X[R==1])
  marmu <- mean(cbind(1,X)%*%lm1$coef)
  
  # Variance estimation
  ipwvar <- diag(VarEst(IPWmrf,ipwpar,data1))[1:2]
  regvar <- diag(VarEst(REGmrf,regpar,data1))[1:2]
  drvar <- diag(VarEst(DRmrf,drpar,data1))[1:2]
  effvar <- diag(VarEst(EFFmrf2,c(effpar,drpar),data1))[1:2]
  
  return(list(mu=mu, ipwmu=ipwpar[1], ipwor=ipwpar[2],
              regmu=regpar[1], regor=regpar[2],
              drmu=drpar[1], dror=drpar[2], 
              effmu=effpar[1], effor=effpar[2], 
              marmu=marmu,
              ipwvar=ipwvar, regvar=regvar, drvar=drvar, effvar=effvar)) 
}

simuTF <- function(N){
  X <- rnorm(N, mean=0,  sd=1) 
  Xr <- cbind(1,X,X^2)
  Xs <- cbind(1,X) 
  Xz <- cbind(1,X,X^2)
  
  dr <- dim(Xr)[2]
  dz <- dim(Xz)[2]
  ds <- dim(Xs)[2]
  
  para = c(-0.5, 0.5, 0.4, 1, 0, 1, -0.3, 0, 1, -0.4, 1, 1)
  gamma <- para[1]
  alpha <- para[2:(ds+1)]
  beta1 <- para[(ds+2):(dr+ds+2)]
  beta2 <- para[(dr+ds+3):(dr+dz+ds+2)]
  beta12 <- beta1[1]
  sigmay2 <- para[dr+dz+ds+3]
  sigmaz2 <- para[dr+dz+ds+4]
  
  # Data-generating process with true propensity and false outcome model
  resp <- plogis(cbind(Xs,Xr,Xz)%*%c(alpha,gamma*beta1[-1],gamma*beta12*beta2) - gamma^2*(sigmay2+beta12^2*sigmaz2)/2)
  R <- rbinom(N,size=1,prob=resp)
  muZ <- Xz%*%beta2 + (R-1)*gamma*beta12*sigmaz2
  Z <- rnorm(N,mean=muZ,sd=sqrt(sigmaz2))
  muY <- cbind(Z,Xr)%*%beta1 + (R-1)*gamma*sigmay2
  Y <- rnorm(N, mean=muY, sd=sqrt(sigmay2))
  mu <- mean(Y)
  Y[R==0] <- NA
  
  data1 <- list(Y=Y,R=R,Z=Z,Xs=cbind(1,X),Xr=cbind(1,X),Xz=cbind(1,X^2))
  
  # IPW estimation 
  ipwpar <- optim(par = c(-0.33,-0.51,0.52,0.41),
                  fn = GMM,
                  f = IPWmrf, data=data1,
                  method = "BFGS", hessian = FALSE)$par
  # REG estimation
  regpar <- optim(par = c(-0.33,-0.51,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                  fn = GMM,
                  f = REGmrf, data=data1,
                  method = "BFGS", hessian = FALSE)$par
  # DR estimation
  drpar <- optim(par = c(-0.33,-0.51,0.52,0.41,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                 fn = GMM,
                 f = DRmrf, data=data1,
                 method = "BFGS", hessian = FALSE)$par
  # EFF estimation
  effres <- EFFmrf(drpar, data1)
  effpar <- apply(effres, 2, mean) + drpar[1:2]
  
  # MAR estimation
  lm1 <- lm(Y[R==1]~X[R==1])
  marmu <- mean(cbind(1,X)%*%lm1$coef)
  
  # Variance estimation
  ipwvar <- diag(VarEst(IPWmrf,ipwpar,data1))[1:2]
  regvar <- diag(VarEst(REGmrf,regpar,data1))[1:2]
  drvar <- diag(VarEst(DRmrf,drpar,data1))[1:2]
  effvar <- diag(VarEst(EFFmrf2,c(effpar,drpar),data1))[1:2]
  
  return(list(mu=mu, ipwmu=ipwpar[1], ipwor=ipwpar[2],
              regmu=regpar[1], regor=regpar[2],
              drmu=drpar[1], dror=drpar[2], 
              effmu=effpar[1], effor=effpar[2], 
              marmu=marmu,
              ipwvar=ipwvar, regvar=regvar, drvar=drvar, effvar=effvar)) 
}


simuFF <- function(N){
  X <- rnorm(N, mean=0,  sd=1) 
  Xr <- cbind(1,X,X^2)
  Xs <- cbind(1,X,X^2) 
  Xz <- cbind(1,X,X^2)
  
  dr <- dim(Xr)[2]
  dz <- dim(Xz)[2]
  ds <- dim(Xs)[2]
  
  para = c(-0.5, 0, 1, 1.5, 1, 0, 1, -0.3, 0, 1, -0.4, 1, 1)
  gamma <- para[1]
  alpha <- para[2:(ds+1)]
  beta1 <- para[(ds+2):(dr+ds+2)]
  beta2 <- para[(dr+ds+3):(dr+dz+ds+2)]
  beta12 <- beta1[1]
  sigmay2 <- para[dr+dz+ds+3]
  sigmaz2 <- para[dr+dz+ds+4]
  
  # Data-generating process with false propensity and false outcome model
  resp <- plogis(cbind(Xs,Xr,Xz)%*%c(alpha,gamma*beta1[-1],gamma*beta12*beta2) - gamma^2*(sigmay2+beta12^2*sigmaz2)/2)
  R <- rbinom(N,size=1,prob=resp)
  muZ <- Xz%*%beta2 + (R-1)*gamma*beta12*sigmaz2
  Z <- rnorm(N,mean=muZ,sd=sqrt(sigmaz2))
  muY <- cbind(Z,Xr)%*%beta1 + (R-1)*gamma*sigmay2
  Y <- rnorm(N, mean=muY, sd=sqrt(sigmay2))
  mu <- mean(Y)
  Y[R==0] <- NA
  
  data1 <- list(Y=Y,R=R,Z=Z,Xs=cbind(1,X),Xr=cbind(1,X),Xz=cbind(1,X^2))
  
  # IPW estimation 
  ipwpar <- optim(par = c(-0.33,-0.51,0.52,0.41),
                  fn = GMM,
                  f = IPWmrf, data=data1,
                  method = "BFGS", hessian = FALSE)$par
  # REG estimation
  regpar <- optim(par = c(-0.33,-0.51,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                  fn = GMM,
                  f = REGmrf, data=data1,
                  method = "BFGS", hessian = FALSE)$par
  # DR estimation
  drpar <- optim(par = c(-0.33,-0.51,0.52,0.41,1.07,0.02,0.98,0.04,-0.47,1.01,1.02),
                 fn = GMM,
                 f = DRmrf, data=data1,
                 method = "BFGS", hessian = FALSE)$par
  # EFF estimation
  effres <- EFFmrf(drpar, data1)
  effpar <- apply(effres, 2, mean) + drpar[1:2]
  
  # MAR estimation
  lm1 <- lm(Y[R==1]~X[R==1])
  marmu <- mean(cbind(1,X)%*%lm1$coef)
  
  # Variance estimation
  ipwvar <- diag(VarEst(IPWmrf,ipwpar,data1))[1:2]
  regvar <- diag(VarEst(REGmrf,regpar,data1))[1:2]
  drvar <- diag(VarEst(DRmrf,drpar,data1))[1:2]
  effvar <- diag(VarEst(EFFmrf2,c(effpar,drpar),data1))[1:2]
  
  return(list(mu=mu, ipwmu=ipwpar[1], ipwor=ipwpar[2],
              regmu=regpar[1], regor=regpar[2],
              drmu=drpar[1], dror=drpar[2], 
              effmu=effpar[1], effor=effpar[2], 
              marmu=marmu,
              ipwvar=ipwvar, regvar=regvar, drvar=drvar, effvar=effvar))  
}




