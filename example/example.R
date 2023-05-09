library(mvtnorm)
library(parallel)
library(mice)
library(miceadds)
library(mitml)
library(lme4)
library(MplusAutomation)
library(survey)
options(survey.lonely.psu = "adjust") # Leite (2016) Chp.2
library(tidyverse)
library(cobalt)



# read-in data ---------------

load("ExampleData.RData")

xnames <- colnames(exdata)[!colnames(exdata)%in%c("trt","y","school","STU_ID" ) ]

ps_x <- paste("trt ~",paste(xnames , collapse=" + "), sep="")
ps_FEx <- paste(ps_x, "+ as.factor(school)")

Rx <- apply(exdata[,xnames],2, function(x){ as.numeric(is.na(x)) })

Rxmiss <- Rx[, colMeans(Rx)>0 ]

Rxnames <- paste0("R.", colnames(Rx[, colMeans(Rx)>0 ]) )
colnames(Rxmiss)<-Rxnames
# cor(Rxmiss)
Rmp <- Rxmiss[ ,c("R.bymathse", "R.bystprep","R.BYSES1")]
Rmpnames <- colnames(Rmp)

ps_FEx_mp <- paste(ps_FEx,"+", paste(Rmpnames , collapse=" + ") )


psmodel <- "ps_FEx" 
weighting <- "clustered"
ymodel <- "yRE" 



# unadjusted
uw.design <- svydesign(ids=~school, # cluster robust standard error
                       nest=TRUE, weights =NULL,
                       data=exdata)
est_unadj <- summary(svyglm(formula = y~trt , family=gaussian(), design=uw.design))$coefficients[2, c(1,2,4)]


# CC -----
dcomp <- exdata[complete.cases(exdata[ ,xnames]), ]


One.CC <- function(dcomp, MIname="CC"){
  
  dimp <- dcomp
  
  # PS model
  psmod_FEx <- glm( formula(ps_FEx), data=dimp, family="binomial")
  psconv_FEx <- (psmod_FEx$converged)
  
  
  psmod <- psmod_FEx
  ps <- fitted(psmod)
  dimp$w <- with(dimp, {
    w <- trt * (1/ps) + (1-trt) * (1/(1-ps)) 
    w
  })
  dat_cl <- with(dimp, {
    wj_sum <- aggregate(data.frame(1, w, w*trt, w*(1-trt)), by=list(school), sum)
    colnames(wj_sum) <- c("school","nj", "wj_sum", "w1j_sum", "w0j_sum")
    merge(dimp, wj_sum, by = "school") -> dat_cl
    
    dat_cl
  })
  wdat_cl <- dat_cl
  wdatcl <- wdat_cl[dat_cl$w0j_sum>0 & dat_cl$w1j_sum>0, ]
  
  wdatcl$scaled_clw <- with(wdatcl, {
    w_sum <- sum(wdatcl$w)
    clw <- (wj_sum/w_sum) * (w*trt/w1j_sum ) + 
      (wj_sum/w_sum) * (w*(1-trt)/w0j_sum )
    scaled_clw <- clw*nrow(wdatcl)/sum(clw) # rescale to facilitate convergence
  })
  
  
  # covariate balance 
  wps <- 1
  cbl.fun <- function(wps){
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights =rep(1, nrow(wdatcl))
    ) -> unblce.clw
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights = wdatcl$scaled_clw
    ) -> blce.clw
    
    clw.cbl <- data.frame(psmod="ps_FEx", weighting="CLW",
                          covname = row.names(blce.clw$Balance),
                          w.smd = blce.clw$Balance$Diff.Adj, u.smd=unblce.clw$Balance$Diff.Adj)
    clw.cbl
  }
  one.cbl.out<-cbl.fun(1)

  # y model 

  datmplus <- with(wdatcl, {
    datmplus <-data.frame(school, trt, y, 
                          scaled_w=scaled_clw)
    datmplus
  })
  
  wfun.Y <- function(datmplus) {
    
    wY.RE <- mplusObject(rdata = datmplus, autov = T, 
                         DATA = "VARIANCES=NOCHECK;",
                         VARIABLE = "
                           USEVARIABLES = trt y scaled_w; 
                           MISSING ARE ALL (999);
                           
                           CLUSTER = school;
                           
                           WITHIN = trt; 
                           BETWEEN = ;
                           WEIGHT = scaled_w; WTSCALE = unscaled;
                           ",
                         ANALYSIS = "
                           Type = Twolevel;
                           Estimator = MLR;
                           ",
                         MODEL = "
                           %WITHIN%
                           y; y on trt;
                           
                           %BETWEEN%
                           y; [y];
                           ",
                         OUTPUT = "NOCHISQUARE TECH3;"
    )
    
    mplusname <- "wY_RE"
    fit.wY.RE <- mplusModeler(wY.RE, dataout = paste0(mplusname, ".dat"),
                              modelout = paste0(mplusname, ".inp"), run = 1) 
    est.wY.RE <- fit.wY.RE$results$parameters$unstandardized[1, c(3, 4, 6)]
    colnames(est.wY.RE) <- c("Est","SE","pval")
    
    clw.est <- data.frame(psmod="ps_FEx", weighting="CLW",
                          data.frame(ymod = c("yRE"), rbind(est.wY.RE))
                          )
    
    clw.est
    
  }
  one.estY.out <- wfun.Y(datmplus)
  
  list(one.cbl.out=one.cbl.out,
       one.estY.out=one.estY.out
  )

}

cc.res.list <- One.CC(dcomp = dcomp, MIname="CC")


# imputation  --------

exdatmi <- exdata %>%
  select( !any_of(c("STU_ID")) )
exdatmi$y <- exdatmi$y/10 # rescale to facilitate convergence


# MI_SL imputation ----
predMat <- make.predictorMatrix(data=exdatmi)
predMat[ ,'school'] <- 0
meth <- make.method(exdatmi) # "pmm"

mice_x <- mice(data = exdatmi, m = 20, method = meth, predictorMatrix = predMat, 
               maxit = 50, seed = 12345 )
plot(mice_x)
complete(mice_x, action = "all", include = FALSE) -> imputedlist.micex


# Each imputed dataset --------
p<-1

imputedlist <- imputedlist.micex
MIname <- "mice_x"

One.Imp <- function(p, imputedlist, MIname){
  
  dimp <- imputedlist[[p]]
  
  # PS model
  psmod_FEx <- glm( formula(ps_FEx), data=dimp, family="binomial")
  psconv_FEx <- (psmod_FEx$converged)
  
  
  dimp_mp <- data.frame(dimp, Rmp)
  
  psmod_FEx_mp <- glm( formula(ps_FEx_mp), data=dimp_mp, family="binomial")
  psconv_FEx_mp <- (psmod_FEx_mp$converged)
  
  
  
  psmod.List <- list(psmod_FEx=psmod_FEx, 
                     psmod_FEx_mp=psmod_FEx_mp
  )
  
  
  psm <- 1
  wfun.psm <- function(psm){
    
    psmod <- psmod.List[[psm]]
    ps <- fitted(psmod)
    
    dimp$w <- with(dimp, {
      w <- trt * (1/ps) + (1-trt) * (1/(1-ps)) 
      w
    })
    
    
    dat_cl <- with(dimp, {
      wj_sum <- aggregate(data.frame(1, w, w*trt, w*(1-trt)), by=list(school), sum)
      colnames(wj_sum) <- c("school","nj", "wj_sum", "w1j_sum", "w0j_sum")
      merge(dimp, wj_sum, by = "school") -> dat_cl
      
      dat_cl
    })
    wdat_cl <- dat_cl
    wdatcl <- wdat_cl[dat_cl$w0j_sum>0 & dat_cl$w1j_sum>0, ]
    
    wdatcl$scaled_clw <- with(wdatcl, {
      w_sum <- sum(wdatcl$w)
      clw <- (wj_sum/w_sum) * (w*trt/w1j_sum ) + 
        (wj_sum/w_sum) * (w*(1-trt)/w0j_sum )
      scaled_clw <- clw * nrow(wdatcl)/sum(clw) 
      scaled_clw
    })
    
    list(wdatcl=wdatcl)
  }
  
  lapply( 1:length(psmod.List), wfun.psm )-> wList
  names(wList) <- names(psmod.List)
  
  
  # covariate balance 
  wps <- 1
  cbl.fun <- function(wps){
    wList[[wps]] -> wdat
    
    psname <- names(wList)[[wps]]
    
    wdatcl <- wdat$wdatcl
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights =rep(1, nrow(wdatcl))
    ) -> unblce.clw
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights = wdatcl$scaled_clw
    ) -> blce.clw
    
    clw.cbl <- data.frame(psmod=psname, weighting="CLW",
                          covname = row.names(blce.clw$Balance),
                          w.smd = blce.clw$Balance$Diff.Adj, u.smd=unblce.clw$Balance$Diff.Adj)
    
    
    rbind(clw.cbl)
  }
  
  lapply( 1:length(psmod.List), cbl.fun) ->cbl.out.list
  
  # y model 
  
  wps <- 1
  estY.fun <- function(wps){
    psname <- names(wList)[[wps]]
    
    wList[[wps]] -> wdat
    
    datmplus <- with(wdat$wdatcl, {
      datmplus <-data.frame(school, trt, y, 
                            scaled_w=scaled_clw)
      datmplus
    })
    
    wfun.Y <- function(datmplus) {
      
      wY.RE <- mplusObject(rdata = datmplus, autov = T, 
                           DATA = "VARIANCES=NOCHECK;",
                           VARIABLE = "
                           USEVARIABLES = trt y scaled_w; 
                           MISSING ARE ALL (999);
                           
                           CLUSTER = school;
                           
                           WITHIN = trt; 
                           BETWEEN = ;
                           WEIGHT = scaled_w; WTSCALE = unscaled;
                           ",
                           ANALYSIS = "
                           Type = Twolevel;
                           Estimator = MLR;
                           ",
                           MODEL = "
                           %WITHIN%
                           y; y on trt;
                           
                           %BETWEEN%
                           y; [y];
                           ",
                           OUTPUT = "NOCHISQUARE TECH3;"
      )
      
      mplusname <- "wY_RE"
      fit.wY.RE <- mplusModeler(wY.RE, dataout = paste0(mplusname, ".dat"),
                                modelout = paste0(mplusname, ".inp"), run = 1) 
      est.wY.RE <- fit.wY.RE$results$parameters$unstandardized[1, c(3, 4, 6)]
      colnames(est.wY.RE) <- c("Est","SE","pval")
      
      data.frame(ymod = c("yRE"), rbind(est.wY.RE))
    }
    
    wfun.Y(datmplus) -> est.clw
    
    clw.est <- data.frame(psmod=psname, weighting="CLW",
                          est.clw)
    
    rbind(clw.est)
    
  }
  
  lapply(1:length(psmod.List), estY.fun) -> estY.out.list
  
  
  cbl.out.df <- do.call(rbind, cbl.out.list) 
  estY.out.df <- do.call(rbind, estY.out.list) 
  
  one.cbl.out <- data.frame(MI=MIname, cbl.out.df) 
  one.estY.out <- data.frame(MI=MIname, estY.out.df) 
  
  list(one.cbl.out=one.cbl.out,
       one.estY.out=one.estY.out
  )
}

micex.res.list <- mclapply(
  1:length(imputedlist.micex), FUN = One.Imp,
  imputedlist = imputedlist.micex, MIname = "mice_x",
  mc.preschedule = T, mc.cores = 20)

# Pool the results  ----
p<-1

# covariate balance
pool.cbl.out <- lapply(1:length(micex.res.list), FUN = function(p){
  data.frame(imp=p, micex.res.list[[p]]$one.cbl.out)
}
)
pool.cbl <- do.call(rbind, pool.cbl.out)

micex.summ.cbl <- pool.cbl %>% 
  group_by(MI, psmod, weighting
  ) %>%
  summarise_at(vars(w.smd, u.smd),
               list( max= ~max(abs(.)), 
                     mean= ~mean(abs(.)) )
  ) %>% relocate(contains("u.smd"), .after = contains("w.smd"))


# estimates ATE
pool.estY.out <- lapply(1:length(micex.res.list), FUN = function(p){
  data.frame(imp=p, micex.res.list[[p]]$one.estY.out)
}
)
pool.estY <- do.call(rbind, pool.estY.out)

micex.estY <- pool.estY %>%
  group_by(MI, psmod, weighting, ymod) %>%
  summarize(
    poolval = testEstimates(qhat = Est, uhat = SE^2 )$estimates[1,c(1,2,5)],
    name = c("Est","SE","pval")
  ) %>%
  pivot_wider(names_from = c("name"), values_from = c("poolval"))



# MI_RE imputation ----

predMat <- make.predictorMatrix(data=exdatmi) 

predMat[ , !colnames(predMat)%in%c("school") ] <- 3 
predMat[ ,'school'] <- -2 # -2 = cluster variable

diag(predMat) <-0

meth <- make.method(exdatmi)
meth[ c("bymathse","byenglse", "byactctl", "byconexp", "byinstmo","bystprep",
        "BYSES1" ) ] <- "2l.pan" 


mice_RIx <- mice( data = exdatmi, m = 20,
                  method = meth, predictorMatrix = predMat,
                  maxit = 50, seed = 12345 )

complete(mice_RIx, action = "all", include = FALSE) -> imputedlist.miceRIx


# Each imputed dataset MI_RE (same as One.Imp) ---------------------
p<-1
imputedlist <- imputedlist.miceRIx
MIname <- "mice_RIx"

miceRIx.res.list <- mclapply(
  1:length(imputedlist.miceRIx), FUN = One.Imp,
  imputedlist = imputedlist.miceRIx, MIname = "mice_RIx",
  mc.preschedule = T, mc.cores = 20)


# Pool the results  -----
p<-1

# covariate balance
pool.cbl.out <- lapply(1:length(miceRIx.res.list), FUN = function(p){
  data.frame(imp=p, miceRIx.res.list[[p]]$one.cbl.out)
}
)
pool.cbl <- do.call(rbind, pool.cbl.out)

miceRIx.summ.cbl <- pool.cbl %>% 
  group_by(MI, psmod, weighting
           # , covname, imp
           ) %>%
  summarise_at(vars(w.smd, u.smd),
               list( max= ~max( abs(.) ), 
                    # n05=~sum(abs(.)>.05), 
                     mean= ~mean( abs(.) ) )
  ) %>% relocate(contains("u.smd"), .after = contains("w.smd"))



# estimates ATE
pool.estY.out <- lapply(1:length(miceRIx.res.list), FUN = function(p){
  data.frame(imp=p, miceRIx.res.list[[p]]$one.estY.out)
}
)
pool.estY <- do.call(rbind, pool.estY.out)

miceRIx.estY <- pool.estY %>%
  group_by(MI, psmod, weighting, ymod) %>%
  summarize(
    poolval = testEstimates(qhat = Est, uhat = SE^2 )$estimates[1,c(1,2,5)],
    name = c("Est","SE","pval")
  ) %>%
  pivot_wider(names_from = c("name"), values_from = c("poolval"))

