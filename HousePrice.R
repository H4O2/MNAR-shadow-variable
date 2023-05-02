##### cfps2010 family house price
###House price of northeast of china
library(stargazer)
library(doParallel)
library(foreach)
source('BasFunctions.R')
load(file = 'cfpsHousePrice.RData')

# Data manipulation
prov1 <- prov%in%c(31,44)
housp[which(housp>1000)] <- 50
Y <- log(housp)
R <- 1 - is.na(Y)
consbpr[which(consbpr<=0)] <- 0.7
Z <- log(consbpr)
flr[flr>0] <- log(flr[flr>0])
X <- cbind(1,prov1,urban,log(dist4), log(siz), log(famsz), flr, log(faminc),recons)
Xr <- Xz <- Xs <- X
data1 <- list(Y=Y,R=R,Z=Z,Xs=Xs,Xr=Xr,Xz=Xz)



### Inverse probability weighting estimator
hmu = 2.6
hgamma = 0.41338589
halpha <- c(1.24650912, -0.59944952, -0.51080139,
            0.12426171, -0.60919234, -0.10941869,  0.21335738,  
            0.23027634, -0.02681936)
ipwpara <- c(hmu, hgamma, halpha)
ipwest <- optim(par = ipwpara,
                fn = GMM,
                f = IPWmrf, data=data1,
                method = "BFGS", hessian = FALSE)$par

### Regression based estimation
# Baseline outcome estimation
lmY <- lm(Y[R==1]~Z[R==1]+X[R==1,]-1)
hbeta1 <- lmY$coef
hsigmay2 <- summary(lmY)$sig^2

# Ancillary outcome estimation
lmZ <- lm(Z[R==1]~X[R==1,]-1)
hbeta2 <- lmZ$coef
hsigmaz2 <- summary(lmZ)$sig^2

regpar <- c(2.59, 0.75, hbeta1, hbeta2, hsigmay2, hsigmaz2)
regest <- optim(par = regpar,
               fn = GMM,
               f = REGmrf, data=data1,
               method = "BFGS", hessian = FALSE)$par

### Doubly robust estimation
drpar <- c(hmu, hgamma, halpha, hbeta1, hbeta2, hsigmay2, hsigmaz2)
drest <- optim(par = drpar,
                fn = GMM,
                f = DRmrf, data=data1,
                method = "BFGS", hessian = FALSE)$par

# Semiparametric efficient estimator
EIF <- EFFmrf(para = drest, data=data1)
effest <- apply(EIF, 2, mean)

### Naive estimators
# Missing at random inverse probability weighting
glmMAR <- glm(R~X-1+Z, family=binomial(link=logit))
ipwMARest <- mean(Y[R==1]/glmMAR$fit[R==1])/length(R)*sum(R)

# Missing at random regression based method
lmMAR <- lm(Y[R==1]~X[R==1,]-1+Z[R==1])
regMAR <- mean(cbind(X,Z)%*%lmMAR$coef)

### Variance estimation
ipwvar <- diag(VarEst(IPWmrf, ipwest, data1))[1:2]
regvar <- diag(VarEst(REGmrf, regest, data1))[1:2]
drvar <- diag(VarEst(DRmrf, drest, data1))[1:2]
effvar <- diag(VarEst2(EFFmrf2, c(effest, drest), data1, EIF = EIF))[1:2]


### Confidence interval
ConfInt(cbind(ipwest[1:2],ipwvar), ci=0.95)
ConfInt(cbind(regest[1:2],regvar), ci=0.95)
ConfInt(cbind(drest[1:2], drvar),  ci=0.95)
ConfInt(cbind(effest,     effvar), ci=0.95)
# Bootstrap estimation of the confidence interval
B = 1000
cl<- makeCluster(detectCores()-2) 
registerDoParallel(cl) 
effBootstrap <- function(data, drpar){
    ind = sample(dim(data$Xr)[1], replace=TRUE)
    Xs2 <- data$Xs[ind,]
    Xr2 <- data$Xr[ind,]
    Xz2 <- data$Xz[ind,]
    Y2 <- data$Y[ind]
    R2 <- data$R[ind]
    Z2 <- data$Z[ind]
    data2 <- list(Y=Y2, R=R2, Z=Z2, Xs=Xs2, Xz=Xz2, Xr=Xr2)
    drest <- optim(par = drpar,
                   fn = GMM,
                   f = DRmrf, data=data2,
                   method = "BFGS", hessian = FALSE)$par
    return(list(effest = drest[1:2] + apply(EFFmrf(para = drest, data=data2), 2, mean)))
}
effestB <- foreach(i=1:B, .inorder=FALSE)%dopar%{
    effBootstrap(data1, drpar)
}
stopCluster(cl)

psi <- gamma <- numeric(B)
for (i in 1:B){
    psi[i] <- effestB[[i]]$effest[1]
    gamma[i] <- effestB[[i]]$effest[2]
}
quantile(psi, c(0.025, 0.975))
quantile(gamma, c(0.025, 0.975))
