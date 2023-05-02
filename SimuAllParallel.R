#Date: 2021/03/01
#### Doubly robust estimation of nonignorable normal data

rm(list=ls())
library(doParallel)
library(foreach)
source('BasFunctions.R')
source('SimuOne.R')

cl<- makeCluster(detectCores()-1) 
registerDoParallel(cl)       
itr <- 1000
gamma <- -0.5

for(mdl in c('FF','FT','TF','TT')){
    for(N in c(500,1500)){
        mu <- ipwmu <- regmu <- drmu <- effmu <- marmu <- numeric(itr)
        ipwor <- regor <- dror <- effor<- numeric(itr)
        ipwvar <- regvar <- drvar <- effvar <- matrix(NA, ncol=2, nrow=itr)
        simuone <- get(paste0('simu', mdl))
         
        simuN <- foreach(i=1:itr, .inorder=FALSE, .packages = "numDeriv")%dopar%{
            simuone(N)
        }
        print(mdl)
        
        for (i in 1:itr){
            mu[i] <- simuN[[i]]$mu
            ipwmu[i] <- simuN[[i]]$ipwmu; 
            regmu[i] <- simuN[[i]]$regmu
            drmu[i] <- simuN[[i]]$drmu
            effmu[i] <- simuN[[i]]$effmu
            marmu[i] <- simuN[[i]]$marmu 
            
            ipwor[i] <- simuN[[i]]$ipwor 
            regor[i] <- simuN[[i]]$regor
            dror[i] <- simuN[[i]]$dror 
            effor[i] <- simuN[[i]]$effor
            
            ipwvar[i,] <- simuN[[i]]$ipwvar
            regvar[i,] <- simuN[[i]]$regvar
            drvar[i,] <- simuN[[i]]$drvar
            effvar[i,] <- simuN[[i]]$effvar
        }
        assign(paste0(mdl,'mu',N), mu)
        assign(paste0(mdl,'ipwmu',N), ipwmu)
        assign(paste0(mdl,'regmu',N), regmu)
        assign(paste0(mdl,'drmu',N), drmu)
        assign(paste0(mdl,'effmu',N), effmu)
        assign(paste0(mdl,'marmu',N), marmu)
        assign(paste0(mdl,'ipwmucvp',N),  mean(na.omit(CoverProb(cbind(ipwmu,ipwvar[,1]), ci=0.95, trvlu=mean(mu)))))     
        assign(paste0(mdl,'regmucvp',N),  mean(na.omit(CoverProb(cbind(regmu,regvar[,1]), ci=0.95, trvlu=mean(mu)))))
        assign(paste0(mdl,'drmucvp',N),   mean(na.omit(CoverProb(cbind(drmu,drvar[,1]),   ci=0.95, trvlu=mean(mu)))))
        assign(paste0(mdl,'effmucvp',N),  mean(na.omit(CoverProb(cbind(effmu,effvar[,1]), ci=0.95, trvlu=mean(mu)))))
        
        assign(paste0(mdl,'ipwor',N), ipwor)
        assign(paste0(mdl,'regor',N), regor)
        assign(paste0(mdl,'dror',N), dror)
        assign(paste0(mdl,'effor',N), effor)
        assign(paste0(mdl,'ipworcvp',N),  mean(na.omit(CoverProb(cbind(ipwor,ipwvar[,2]), ci=0.95, trvlu=gamma))))     
        assign(paste0(mdl,'regorcvp',N),  mean(na.omit(CoverProb(cbind(regor,regvar[,2]), ci=0.95, trvlu=gamma))))     
        assign(paste0(mdl,'drorcvp',N),   mean(na.omit(CoverProb(cbind(dror,drvar[,2]),   ci=0.95, trvlu=gamma))))  
        assign(paste0(mdl,'efforcvp',N),  mean(na.omit(CoverProb(cbind(effor,effvar[,2]), ci=0.95, trvlu=gamma))))    
    }
}

stopCluster(cl)
save.image(file='SimuResults.RData')

# summary of simulations
load("Simu500.RData")
load("Simu1000.RData")
for(mdl in  c('FF','FT','TF','TT')){
    
    # boxplot mu
    pdf(paste0('simu',mdl,'mu.pdf'), width=10, height=6.3)
    
    boxplot(get(paste0(mdl,'effmu500')),get(paste0(mdl,'effmu1000')),
            get(paste0(mdl,'drmu500')),get(paste0(mdl,'drmu1000')),
            get(paste0(mdl,'ipwmu500')),get(paste0(mdl,'ipwmu1000')),
            get(paste0(mdl,'regmu500')),get(paste0(mdl,'regmu1000')),
            get(paste0(mdl,'marmu500')),get(paste0(mdl,'marmu1000')),
            col=c('gray93','gray','gray93','gray',
                  'gray93','gray','gray93','gray','gray93','gray'),
            outline=FALSE, 
            xaxt='n', 
            at=c(1,2, 3.5, 4.5, 6,7, 8.5, 9.5,11,12), 
            boxwex=0.5, cex.main=1.5, lwd=1.25,
            cex.axis=2)
    axis(1, at=seq(1.5,by=2.5,length.out=5), 
         labels=c('EFF', 'DR', 'IPW', 'REG','MAR'), cex.axis=2)
    abline(h=mean(get(paste0(mdl,'mu1000'))), col='gray60')
    dev.off()
    
    # boxplot mu
    pdf(paste0('simu',mdl,'mu2.pdf'), width=10, height=6.3)
    
    boxplot(get(paste0(mdl,'effmu500')),get(paste0(mdl,'effmu1000')),
            get(paste0(mdl,'drmu500')),get(paste0(mdl,'drmu1000')),
            get(paste0(mdl,'ipwmu500')),get(paste0(mdl,'ipwmu1000')),
            get(paste0(mdl,'regmu500')),get(paste0(mdl,'regmu1000')),
            col=c('gray93','gray','gray93','gray',
                  'gray93','gray','gray93','gray'),
            outline=FALSE, 
            xaxt='n', 
            at=c(1,2, 3.5, 4.5, 6,7, 8.5, 9.5), 
            boxwex=0.5, cex.main=1.5, lwd=1.25,
            cex.axis=2)
    axis(1, at=seq(1.5,by=2.5,length.out=4), 
         labels=c('EFF', 'DR', 'IPW', 'REG'), cex.axis=2)
    abline(h=mean(get(paste0(mdl,'mu1000'))), col='gray60')
    dev.off()
    
    
    # boxplot or
    pdf(paste0('simu',mdl,'or.pdf'),width=10,height=6.3)
    
    boxplot(get(paste0(mdl,'effor500')),get(paste0(mdl,'effor1000')),
            get(paste0(mdl,'dror500')),get(paste0(mdl,'dror1000')),
            get(paste0(mdl,'ipwor500')),get(paste0(mdl,'ipwor1000')),
            get(paste0(mdl,'regor500')),get(paste0(mdl,'regor1000')),
            col=c('gray93','gray','gray93','gray',
                  'gray93','gray','gray93','gray'),
            outline=FALSE, 
            xaxt='n',
            at=c(1,2, 3.5, 4.5, 6,7, 8.5, 9.5),
            boxwex=0.55, cex.main=2, lwd=1.25,
            cex.axis=2)
    axis(1, at=seq(1.5,by=2.5,length.out=4), 
         labels=c('EFF','DR','IPW','REG'), cex.axis=2)
    abline(h=gamma, col='gray60')
    dev.off()
    
    # coverage
    for(N in c(500,1000)){
        assign(paste0(mdl,'cvp',N), c(get(paste0(mdl,'effmucvp',N)), get(paste0(mdl,'drmucvp',N)),get(paste0(mdl,'ipwmucvp',N)),get(paste0(mdl,'regmucvp',N)), get(paste0(mdl,'efforcvp',N)), get(paste0(mdl,'drorcvp',N)),get(paste0(mdl,'ipworcvp',N)),get(paste0(mdl,'regorcvp',N))))
    }
}

allcvp <- rbind(FTcvp500,FTcvp1000,TFcvp500,TFcvp1000,
                TTcvp500,TTcvp1000,FFcvp500,FFcvp1000)
dimnames(allcvp) <- list(c('FT','FT','TF','TF','TT','TT','FF','FF'),
                         c('EFF','DR','IPW','REG','EFF','DR','IPW','REG'))
write.table(allcvp,'cvp.txt', sep='\t')
 
