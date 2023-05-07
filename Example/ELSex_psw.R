library(mvtnorm)
library(parallel)
library(mice)
library(miceadds)
library(mitml)
library(lme4)
library(MplusAutomation)
library(survey)

library(tidyverse)
library(cobalt)
# survey package for cluster-robust SE with single-level outcome model # following Leite (2017)
#source("~/Library/CloudStorage/Box-Box/Labs/PS_Missingdata/Simulations/fun_r1methods_sandv.R")

#load("~/Library/CloudStorage/Box-Box/Labs/PS_Missingdata/Example/ELS_career_academy.RData")
#write.csv(ELS.data, "ELS.data.csv", row.names = F)

# Ignore sampling weights for illustration (as in Leite et al., 2021, missing data PS) 
# surveyDesign1 <- svydesign(ids=~psu, strata=~STRAT_ID, weights=~bystuwt, data = exELS, nest=T)

# Leite et al. (2021) example
# From the 31 covariates included in analysis, only 12 covariates were selected because they had at least 5% of missing data. This was done in this demonstration to increase the overall percentage of missing data to 16.9%.


options(survey.lonely.psu = "adjust") # Leite (2016) Chp.2
# help(svyCprod) 


# read-in data ---------------
load("exdatdum.RData")

# exdatdum %>% group_by(school) %>%
#   summarise_all( ~var(., na.rm=T) )->clusvar
# covariateNames.L2 <- colnames(exdatdum)[colSums(clusvar, na.rm = T)==0]

covariateNames.L2 <- c("byacclim", "bysctrl",  "byurban",  "byregion", "BY10FLP" )

xnames <- colnames(exdatdum)[
  !colnames(exdatdum)%in%c("trt","y","school","STU_ID",
                           covariateNames.L2) ]

ps_x <- paste("trt ~",paste(xnames , collapse=" + "), sep="")
ps_FEx <- paste(ps_x, "+ as.factor(school)")
# no L2 covariates, so the complete-cases for ps_FEx and ps_REx are the same
ps_RIx <- paste(ps_x, # paste(covariateNames.L2, collapse="+"),
                "+ ( 1 | school )")

cmxNames <- paste0("cm.",xnames)
ps_RIxcmx <- paste(ps_RIx,"+", paste(cmxNames , collapse=" + ") )

Rx <- apply(exdatdum[,xnames],2, function(x){ as.numeric(is.na(x)) })

summary(colMeans(Rx)[colMeans(Rx)>0])

Rxmiss <- Rx[, colMeans(Rx)>0 ]

Rxnames <- paste0("R.", colnames(Rx[, colMeans(Rx)>0 ]) )
colnames(Rxmiss)<-Rxnames
apply(cor(Rxmiss), 2, function(x){ max(abs(x)[abs(x)!=1]) })
Rmp <- Rxmiss[ ,c("R.bymathse", "R.bystprep","R.BYSES1")]
Rmpnames <- colnames(Rmp)
ps_RIx_mp <- paste(ps_RIx,"+", paste(Rmpnames , collapse=" + ") )

ps_FEx_mp <- paste(ps_FEx,"+", paste(Rmpnames , collapse=" + ") )


psmodel <- "ps_FEx" # psmodel <- "ps_RIx"
weighting <- "clustered"
ymodel <- "ySl" # ymodel<- "yRI"



# unadjusted ----

uw.design <- svydesign(ids=~school, nest=TRUE, weights =NULL,
                       data=exdatdum)
est_unadj <- summary( svyglm(formula = y~trt , family=gaussian(), design=uw.design))$coefficients[2, c(1,2,4)]
# summary(lm( y~trt, data = exdatdum ))


# adjusted -----
exIPW.Outcome <- function( imps ){
  with(imps, {
    
    if( psmodel=="ps_FEx" ){ 
      psmod <- glm( formula(ps_FEx), family="binomial") 
      psconv <- (psmod$converged)
    }
    
    if( psmodel=="ps_RIx" ){ 
      psmod <- glmer( ps_RIx, family="binomial"  ) 
      psconv <- !any( grepl("failed to converge", psmod@optinfo$conv$lme4$messages) )
    }
    
    ps <- fitted(psmod)
    w <- trt * (1/ps) + (1-trt) * (1/(1-ps))
    
    if( weighting=="IPW" ){
      scaled_w <- (trt==1)*w /sum( w[trt ==1]) +  (trt==0)*w /sum( w[trt ==0])
    }
    
    if( weighting=="clustered" ){
      wj_sum <- aggregate(data.frame(1, w, w*trt, w*(1-trt)), by=list(school), sum)
      colnames(wj_sum) <- c("school","nj", "wj_sum", "w1j_sum", "w0j_sum")
      merge(data.frame(y=y, trt=trt, school=school, scaled_w=w ),wj_sum, by = "school") -> dat_cl
      datcl <- dat_cl[dat_cl$w0j_sum>0 & dat_cl$w1j_sum>0, ]
      w_sum <-sum( wj_sum$wj_sum[ wj_sum$school %in% datcl$school ] )
      
      datcl$cl_w <-  ifelse(
        datcl$trt==1, 
        (datcl$wj_sum/w_sum) *(datcl$scaled_w* datcl$trt / datcl$w1j_sum ),
        (datcl$wj_sum/w_sum) *(datcl$scaled_w* (1-datcl$trt) / datcl$w0j_sum )
      )
      
      
    }
    
    if( weighting!="clustered"){
      datmplus <-data.frame( school, trt, y, scaled_w)
      datmplus$scaled_w <- datmplus$scaled_w * nrow(datmplus)/sum(datmplus$scaled_w)
      psdf <- model.frame(psmod)
      zdf <- data.frame( bysctrl, byurban, byregion)
    }
    
    if( weighting=="clustered" ){
      datmplus <-data.frame( school=datcl$school, trt=datcl$trt, y=datcl$y, scaled_w=datcl$cl_w)
      datmplus$scaled_w <- datmplus$scaled_w * nrow(datmplus)/sum(datmplus$scaled_w)
      dd <- model.frame(psmod)
      psdf <- dd[ dat_cl$w0j_sum>0 & dat_cl$w1j_sum>0,  ]
      zdf <- data.frame(bysctrl, byurban, byregion)[ dat_cl$w0j_sum>0 & dat_cl$w1j_sum>0,  ]
    }
    
    if( ymodel=="yRI" ){
      
      
      pid <- Sys.getpid()
      write.table(datmplus, paste0(pid, "datmplus.txt"), 
                  na="999", quote = F,row.names = F,col.names = F )
      
      cat(
        sep = " ",
        " Data: 
        file =", paste0(pid, "datmplus.txt" ) ,";
        VARIANCES=NOCHECK;
        
        Variable:
        NAMES = ", colnames(datmplus), ";
        USEVARIABLES = trt y scaled_w; 
        MISSING ARE ALL (999);
        CLUSTER = school;
        
        WITHIN = trt; !not model trt to be latent cluster mean centered; 
        BETWEEN = ;
        WEIGHT = scaled_w; WTSCALE = unscaled;
        
        !Define: trt (groupmean); !not groupmean center t;
        
        Analysis:
        Type = Twolevel;
        Estimator = MLR;
        
        Model:
        %WITHIN%
        y; y on trt;
        
        %BETWEEN%
        y; [y];
        
        Output:
        NOCHISQUARE TECH1 TECH3;
        ", 
        file = paste0( "wY_RI.inp")
      )
      runModels( paste0( "wY_RI.inp" ), logFile = NULL)
      wY_RI_out <- readModels( paste0(  "wY_RI.out"), what = c("all"))
      est <- wY_RI_out$parameters$unstandardized[1, c(3, 4, 6)]
    }
    
    if( ymodel == "ySl"){
      
      psw.design <- svydesign(ids=~school, nest=TRUE , weights =~scaled_w,data=datmplus)
      est <-summary( svyglm(formula = y~trt, family=gaussian(), design=psw.design))$coefficients[2, c(1,2,4)]
      
    }
    
    out <- data.frame(ps=psmodel, ymodel=ymodel, weighting=weighting
                      ,matrix(est,nrow = 1), psconv=psconv) 
    colnames(out) <- c("ps","ymodel","weighting","Est","SE","pval","psconv")
    
    outall <- mget( ls(envir = environment()) )
    outall
  }) -> outall
  
  outall
}


w.Y.RI <- function(datmplus){
  
  write.table(datmplus, paste0("datmplus.txt"), 
              na="999", quote = F,row.names = F,col.names = F )
  
  cat(
    sep = " ",
    " Data: 
    file =", paste0( "datmplus.txt" ) ,";
    VARIANCES=NOCHECK;
    
    Variable:
    NAMES = ", colnames(datmplus), ";
    USEVARIABLES = trt y scaled_w; 
    MISSING ARE ALL (999);
    CLUSTER = school;
    
    WITHIN = trt; !not model trt to be latent cluster mean centered; 
    BETWEEN = ;
    WEIGHT = scaled_w; WTSCALE = unscaled;
    
    !Define: trt (groupmean); !not groupmean center t;
    
    Analysis:
    Type = Twolevel;
    Estimator = MLR;
    
    Model:
    %WITHIN%
    y; y on trt;
    
    %BETWEEN%
    y; [y];
    
    Output:
    NOCHISQUARE TECH1 TECH3;
    ", 
    file = paste0( "wY_RI.inp")
  )
  runModels( paste0( "wY_RI.inp" ), logFile = NULL)
  wY_RI_out <- readModels( paste0(  "wY_RI.out"), what = c("all"))
  est <- wY_RI_out$parameters$unstandardized[1, c(3, 4, 6)]
  
  est
}

library(cobalt)

# CC -----
dcomp <- exdatdum[complete.cases(exdatdum[ ,xnames]), ]


dcomp %>% group_by(school) %>%
  mutate( nj1 = sum(trt),
          nj0 = sum(1-trt),
          nj = n()) -> dcomp
dcomp %>% group_by(school) %>%
  summarise( nj1 = sum(trt),
             nj0 = sum(1-trt),
             nj = n()) -> dcompsum

summary(dcompsum$nj)


#  CC dataset --------
One.CC <- function(dcomp, MIname="CC"){
  
  dimp <- dcomp
  
  # PS model
  psmod_FEx <- glm( formula(ps_FEx), data=dimp, family="binomial")
  psconv_FEx <- (psmod_FEx$converged)
  
  psmod_RIx <- glmer( ps_RIx, data=dimp,family="binomial"  ) 
  psconv_RIx <- !any( grepl("failed to converge", 
                            psmod_RIx@optinfo$conv$lme4$messages) )
  
  
  cmx <- aggregate(dimp[ ,xnames], by = list(dimp$school), FUN = mean, na.rm=T)
  colnames(cmx) <- c("school", cmxNames)
  cor(cmx[,-1])
  dimp_cmx <- merge(dimp, cmx, by = "school")
  
  psmod_RIxcmx <- glmer( ps_RIxcmx, data=dimp_cmx,family="binomial"  ) 
  psconv_RIxcmx <- !any( grepl("failed to converge", 
                               psmod_RIxcmx@optinfo$conv$lme4$messages) )
  
  
  # weighting
  psmod.List <- list(psmod_FEx=psmod_FEx, psmod_RIx=psmod_RIx, 
                     psmod_RIxcmx=psmod_RIxcmx
  )
  
  psconv.List <- list(psconv_FEx=psconv_FEx, psconv_RIx=psconv_RIx, 
                      psconv_RIxcmx=psconv_RIxcmx
  )
  
  psm <- 1
  wfun.psm <- function(psm){
    
    psmod <- psmod.List[[psm]]
    ps <- fitted(psmod)
    
    dimp$w <- with(dimp, {
      w <- trt * (1/ps) + (1-trt) * (1/(1-ps)) 
      w
    })
    dimp$scaled_ipw <- with(dimp, { 
      scaled_w <- (trt==1)*w /sum( w[trt ==1]) +  (trt==0)*w /sum( w[trt ==0])
      scaled_ipw <- scaled_w * nrow(dimp)/sum(scaled_w)
      scaled_ipw
    })
    
    # which(exdatdum[,"y"] != dimp$y)
    wdimp <- dimp
    
    dat_cl <- with(dimp, {
      wj_sum <- aggregate(data.frame(1, w, w*trt, w*(1-trt)), by=list(school), sum)
      colnames(wj_sum) <- c("school","nj", "wj_sum", "w1j_sum", "w0j_sum")
      merge(dimp, wj_sum, by = "school") -> dat_cl
      
      dat_cl
    })
    # which(dat_cl[,"y"] != exdatdum[,"y"])
    wdat_cl <- dat_cl
    wdatcl <- wdat_cl[dat_cl$w0j_sum>0 & dat_cl$w1j_sum>0, ]
    
    wdatcl$scaled_clw <- with(wdatcl, {
      w_sum <- sum(wdatcl$w)
      clw <- (wj_sum/w_sum) * (w*trt/w1j_sum ) + 
        (wj_sum/w_sum) * (w*(1-trt)/w0j_sum )
      
      scaled_clw <- clw * nrow(wdatcl)/sum(clw)
      scaled_clw
    })
    
    
    list(wdimp=wdimp, wdatcl=wdatcl)
  }
  
  lapply( 1:length(psmod.List), wfun.psm )-> wList
  names(wList) <- names(psmod.List)
  
  
  # covariate balance 
  wps <- 1
  cbl.fun <- function(wps){
    wList[[wps]] -> wdat
    
    psname <- names(wList)[[wps]]
    
    znames <- c("bysctrl",  "byurban",  "byregion")
    
    # ipw
    wdimp <- wdat$wdimp
    
    bal.tab(wdimp$trt ~ wdimp[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights =rep(1, nrow(wdimp))
    ) -> unblce.ipw
    
    bal.tab(wdimp$trt ~ wdimp[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights = wdimp$scaled_ipw
    ) -> blce.ipw
    
    ipw.cbl <- data.frame(psmod=psname, weighting="IPW",
                          covname = row.names(blce.ipw$Balance),
                          w.smd = blce.ipw$Balance$Diff.Adj, u.smd=unblce.ipw$Balance$Diff.Adj)
    
    
    # clw
    wdatcl <- wdat$wdatcl
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights =rep(1, nrow(wdatcl))
    ) -> unblce.clw
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights = wdatcl$scaled_clw
    ) -> blce.clw
    
    clw.cbl <- data.frame(psmod=psname, weighting="clustered",
                          covname = row.names(blce.clw$Balance),
                          w.smd = blce.clw$Balance$Diff.Adj, u.smd=unblce.clw$Balance$Diff.Adj)
    
    
    rbind(clw.cbl, ipw.cbl)
  }
  
  #cbl.fun(1) ->cbl.out
  lapply( 1:length(psmod.List), cbl.fun) ->cbl.out.list
  
  # RI for a covariate
  estPre.fun <- function(wps){
    psname <- names(wList)[[wps]]
    
    wList[[wps]] -> wdat
    
    predatmplus.ipw <- with(wdat$wdimp, {
      datmplus <-data.frame(school, trt, pre=bytxmstd, 
                            scaled_w=scaled_ipw)
      datmplus
    })
    
    predatmplus.clw <- with(wdat$wdatcl, {
      datmplus <-data.frame(school, trt, pre=bytxmstd, 
                            scaled_w=scaled_clw)
      datmplus
    })
    
    wfun.pre <- function(predatmplus) {
      
      wpre.RI <- mplusObject(rdata = predatmplus, autov = T, 
                             DATA = "VARIANCES=NOCHECK;",
                             VARIABLE = "
 USEVARIABLES = trt pre scaled_w; 
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
 pre; pre on trt;
 
 %BETWEEN%
 pre; [pre];
 ",
                             OUTPUT = "NOCHISQUARE TECH3;"
      )
      
      mplusname <- tempfile("wpre_RI", tmpdir = "tmp")
      fit.wpre.RI <- mplusModeler(wpre.RI, 
                                  dataout = paste0(mplusname, ".dat"),
                                  modelout = paste0(mplusname, ".inp"),
                                  run = 1) 
      
      est.wpre.RI <- fit.wpre.RI$results$parameters$unstandardized[1, c(3, 4, 6)]
      colnames(est.wpre.RI) <- c("Est","SE","pval")
      
      # pre_SL
      prepsw.design <- svydesign(ids=~school, nest=TRUE, 
                              weights =~scaled_w, 
                              data=predatmplus)
      est.wpre.SL <-summary(svyglm(formula = pre ~ trt, family=gaussian(),
                                   design=prepsw.design))$coefficients[2, c(1,2,4)]
      names(est.wpre.SL) <- c("Est","SE","pval")
      
      data.frame(premod = c("preRI","preSL"), rbind(est.wpre.RI, est.wpre.SL))
    }
    
    
    wfun.pre(predatmplus.ipw) -> preest.ipw
    wfun.pre(predatmplus.clw) -> preest.clw
    
    ipw.preest <- data.frame(psmod=psname, weighting="IPW",
                             preest.ipw)
    
    clw.preest <- data.frame(psmod=psname, weighting="clustered",
                             preest.clw)
    
    rbind(clw.preest, ipw.preest)
    
  }
  lapply(1:length(psmod.List), estPre.fun) -> estpre.out.list
  
  
  # y model 
  
  wps <- 1
  estY.fun <- function(wps){
    psname <- names(wList)[[wps]]
    
    wList[[wps]] -> wdat
    
    datmplus.ipw <- with(wdat$wdimp, {
      datmplus <-data.frame(school, trt, y, 
                            scaled_w=scaled_ipw)
      datmplus
    })
    
    datmplus.clw <- with(wdat$wdatcl, {
      datmplus <-data.frame(school, trt, y, 
                            scaled_w=scaled_clw)
      datmplus
    })
    
    wfun.Y <- function(datmplus) {
      
      wY.RI <- mplusObject(rdata = datmplus, autov = T, 
                           DATA = "VARIANCES=NOCHECK;",
                           VARIABLE = "
                           USEVARIABLES = trt y scaled_w; 
                           MISSING ARE ALL (999);
                           
                           CLUSTER = school;
                           
                           WITHIN = trt; !not model trt to be latent cluster mean centered; 
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
      
      mplusname <- tempfile("wY_RI", tmpdir = "tmp")
      fit.wY.RI <- mplusModeler(wY.RI, dataout = paste0(mplusname, ".dat"),
                                modelout = paste0(mplusname, ".inp"), run = 1) # If "run" greater than one, the model is bootstrapped with run replications as well as the basic model.
      #fit.wY.RI <- readModels( paste0(  "wY_RI.out"), what = c("all"))
      
      est.wY.RI <- fit.wY.RI$results$parameters$unstandardized[1, c(3, 4, 6)]
      colnames(est.wY.RI) <- c("Est","SE","pval")
      
      # Y_SL
      psw.design <- svydesign(ids=~school, nest=TRUE, 
                              weights =~scaled_w, data=datmplus)
      est.wY.SL <-summary(svyglm(formula = y~trt, family=gaussian(), design=psw.design))$coefficients[2, c(1,2,4)]
      names(est.wY.SL) <- c("Est","SE","pval")
      
      data.frame(ymod = c("yRI","ySL"), rbind(est.wY.RI, est.wY.SL))
    }
    
    
    wfun.Y(datmplus.ipw) -> est.ipw
    wfun.Y(datmplus.clw) -> est.clw
    
    ipw.est <- data.frame(psmod=psname, weighting="IPW",
                          est.ipw)
    
    clw.est <- data.frame(psmod=psname, weighting="clustered",
                          est.clw)
    
    rbind(clw.est, ipw.est)
    
  }
  
  lapply(1:length(psmod.List), estY.fun) -> estY.out.list
  
  
  cbl.out.df <- do.call(rbind, cbl.out.list) 
  estpre.out.df <- do.call(rbind, estpre.out.list) 
  estY.out.df <- do.call(rbind, estY.out.list) 
  
  one.cbl.out <- data.frame(MI=MIname, cbl.out.df) 
  one.estpre.out <- data.frame(MI=MIname, estpre.out.df) 
  one.estY.out <- data.frame(MI=MIname, estY.out.df) 
  
  list(one.cbl.out=one.cbl.out,
       one.estpre.out=one.estpre.out, 
       one.estY.out=one.estY.out,
       psconv.List=psconv.List
  )
}

cc.res.list <- One.CC(dcomp = dcomp, MIname="CC")

save(cc.res.list, file = "cc.res.list.RData")

load("~/Library/CloudStorage/Box-Box/Labs/PS_Missingdata/Example/cc.res.list.RData")


ccx.res.list <- list(cc.res.list)
p<-1

# covariate balance
pool.cbl.out <- lapply(1:length(ccx.res.list), FUN = function(p){
  data.frame(imp=p, ccx.res.list[[p]]$one.cbl.out)
}
)
pool.cbl <- do.call(rbind, pool.cbl.out)

ccx.summ.cbl <- pool.cbl %>% 
  group_by(MI, psmod, weighting
           # , covname, imp
  ) %>%
  summarise_at(vars(w.smd, u.smd),
               list( max= ~max(abs(.)), 
                     # n05=~sum(abs(.)>.05), 
                     mean= ~mean(abs(.)) )
  ) %>% relocate(contains("u.smd"), .after = contains("w.smd"))



# BY level-2 covariates
ccx.summ.cbl.L2L1 <- pool.cbl %>% 
  mutate( covLevel = ifelse(covname %in% c("bysctrl",  "byurban",  "byregion"), "L2", "L1")  ) %>%
  group_by(MI, psmod, weighting
           ,covLevel # , covname
  ) %>%
  summarise_at(vars(w.smd, u.smd),
               list( max= ~max(abs(.)), 
                     # n05=~sum(abs(.)>.05), 
                     mean= ~mean(abs(.)) )
  ) %>% relocate(contains("u.smd"), .after = contains("w.smd"))



# estimates ATE
pool.estY.out <- lapply(1:length(ccx.res.list), FUN = function(p){
  data.frame(imp=p, ccx.res.list[[p]]$one.estY.out)
}
)
pool.estY <- do.call(rbind, pool.estY.out)

ccx.estY <- pool.estY[,-1]



# # estimates pre
# pool.estpre.out <- lapply(1:length(ccx.res.list), FUN = function(p){
#   data.frame(imp=p, ccx.res.list[[p]]$one.estpre.out)
# }
# )
# pool.estpre <- do.call(rbind, pool.estpre.out)
# 
# ccx.estpre <- pool.estpre %>%
#   group_by(MI, psmod, weighting, premod) %>%
#   summarize(
#     poolval = testEstimates(qhat = Est, uhat = SE^2 )$estimates[1,c(1,2,5)],
#     name = c("Est","SE","pval")
#             ) %>%
#   pivot_wider(names_from = c("name"), values_from = c("poolval"))


# imputation prep data --------
# make dummies "Note that categorical variables must be converted into a numeric data type beforehand in order to be imputed using mice." (Grund et al., 2018)
#model.matrix(formula(ps_x), data=exdat) -> xdat

exdatmi <- exdatdum %>%
  select( !any_of(c(covariateNames.L2, "STU_ID")) )

colnames(exdatmi) <- gsub(">=","", colnames(exdatmi))
var(exdatmi$y)
apply(exdatmi, 2, function(x){var(x, na.rm = T)})

save(exdatmi,file = "exdatmi.RData")

load("exdatmi.RData")


# MI_SL imputation ----
predMat <- make.predictorMatrix(data=exdatmi)
predMat[ ,'school'] <- 0
meth <- make.method(exdatmi) # "pmm"

mice_x <- mice(data = exdatmi, m = 20, method = meth, predictorMatrix = predMat, 
               maxit = 50, seed = 12345 )

save(mice_x, file = "Eg_mice_x.RData")
plot(mice_x)

mice_x$method
mice_x$predictorMatrix

load("Eg_mice_x.RData")
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
  
  psmod_RIx <- glmer( ps_RIx, data=dimp,family="binomial"  ) 
  psconv_RIx <- !any( grepl("failed to converge", 
                            psmod_RIx@optinfo$conv$lme4$messages) )
  
  dimp_mp <- data.frame(dimp, Rmp)
  
  psmod_FEx_mp <- glm( formula(ps_FEx_mp), data=dimp_mp, family="binomial")
  psconv_FEx_mp <- (psmod_FEx_mp$converged)
  
  psmod_RIx_mp <- glmer( ps_RIx_mp, data=dimp_mp,family="binomial"  ) 
  psconv_RIx_mp <- !any( grepl("failed to converge", 
                               psmod_RIx_mp@optinfo$conv$lme4$messages) )
  
  
  cmx <- aggregate(dimp[ ,xnames], by = list(dimp$school), FUN = mean, na.rm=T)
  colnames(cmx) <- c("school", cmxNames)
  cor(cmx[,-1])
  dimp_cmx <- merge(dimp, cmx, by = "school")
  
  psmod_RIxcmx <- glmer( ps_RIxcmx, data=dimp_cmx,family="binomial"  ) 
  psconv_RIxcmx <- !any( grepl("failed to converge", 
                               psmod_RIxcmx@optinfo$conv$lme4$messages) )
  
  
  # weighting
  psmod.List <- list(psmod_FEx=psmod_FEx, psmod_RIx=psmod_RIx, 
                     psmod_FEx_mp=psmod_FEx_mp, psmod_RIx_mp=psmod_RIx_mp,
                     psmod_RIxcmx=psmod_RIxcmx
  )
  
  psconv.List <- list(psconv_FEx=psconv_FEx, psconv_RIx=psconv_RIx, 
                      psconv_FEx_mp=psconv_FEx_mp, psconv_RIx_mp=psconv_RIx_mp,
                      psconv_RIxcmx=psconv_RIxcmx
  )
  
  psm <- 1
  wfun.psm <- function(psm){
    
    psmod <- psmod.List[[psm]]
    ps <- fitted(psmod)
    
    dimp$w <- with(dimp, {
      w <- trt * (1/ps) + (1-trt) * (1/(1-ps)) 
      w
    })
    dimp$scaled_ipw <- with(dimp, { 
      scaled_w <- (trt==1)*w /sum( w[trt ==1]) +  (trt==0)*w /sum( w[trt ==0])
      scaled_ipw <- scaled_w * nrow(dimp)/sum(scaled_w)
      scaled_ipw
    })
    
    # which(exdatdum[,"y"] != dimp$y)
    wdimp <- data.frame(dimp, exdatdum[,c("bysctrl",  "byurban",  "byregion")])
    
    dat_cl <- with(dimp, {
      wj_sum <- aggregate(data.frame(1, w, w*trt, w*(1-trt)), by=list(school), sum)
      colnames(wj_sum) <- c("school","nj", "wj_sum", "w1j_sum", "w0j_sum")
      merge(dimp, wj_sum, by = "school") -> dat_cl
      
      dat_cl
    })
    # which(dat_cl[,"y"] != exdatdum[,"y"])
    wdat_cl <- data.frame(dat_cl, exdatdum[,c("bysctrl",  "byurban",  "byregion")])
    wdatcl <- wdat_cl[dat_cl$w0j_sum>0 & dat_cl$w1j_sum>0, ]
    
    wdatcl$scaled_clw <- with(wdatcl, {
      w_sum <- sum(wdatcl$w)
      clw <- (wj_sum/w_sum) * (w*trt/w1j_sum ) + 
        (wj_sum/w_sum) * (w*(1-trt)/w0j_sum )
      
      scaled_clw <- clw * nrow(wdatcl)/sum(clw)
      scaled_clw
    })
    
    
    list(wdimp=wdimp, wdatcl=wdatcl)
  }
  
  lapply( 1:length(psmod.List), wfun.psm )-> wList
  names(wList) <- names(psmod.List)
  
  
  # covariate balance 
  wps <- 1
  cbl.fun <- function(wps){
    wList[[wps]] -> wdat
    
    psname <- names(wList)[[wps]]
    
    znames <- c("bysctrl",  "byurban",  "byregion")
    
    # ipw
    wdimp <- wdat$wdimp
    
    bal.tab(wdimp$trt ~ wdimp[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights =rep(1, nrow(wdimp))
    ) -> unblce.ipw
    
    bal.tab(wdimp$trt ~ wdimp[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights = wdimp$scaled_ipw
    ) -> blce.ipw
    
    ipw.cbl <- data.frame(psmod=psname, weighting="IPW",
                          covname = row.names(blce.ipw$Balance),
                          w.smd = blce.ipw$Balance$Diff.Adj, u.smd=unblce.ipw$Balance$Diff.Adj)
    
    
    # clw
    wdatcl <- wdat$wdatcl
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights =rep(1, nrow(wdatcl))
    ) -> unblce.clw
    
    bal.tab(wdatcl$trt ~ wdatcl[ ,c(xnames, znames)],
            s.d.denom = "pooled", continuous="std", binary="std",
            weights = wdatcl$scaled_clw
    ) -> blce.clw
    
    clw.cbl <- data.frame(psmod=psname, weighting="clustered",
                          covname = row.names(blce.clw$Balance),
                          w.smd = blce.clw$Balance$Diff.Adj, u.smd=unblce.clw$Balance$Diff.Adj)
    
    
    rbind(clw.cbl, ipw.cbl)
  }
  
  #cbl.fun(1) ->cbl.out
  lapply( 1:length(psmod.List), cbl.fun) ->cbl.out.list
  
  # RI for a covariate
  estPre.fun <- function(wps){
    psname <- names(wList)[[wps]]
    
    wList[[wps]] -> wdat
    
    predatmplus.ipw <- with(wdat$wdimp, {
      datmplus <-data.frame(school, trt, pre=bytxmstd, 
                            scaled_w=scaled_ipw)
      datmplus
    })
    
    predatmplus.clw <- with(wdat$wdatcl, {
      datmplus <-data.frame(school, trt, pre=bytxmstd, 
                            scaled_w=scaled_clw)
      datmplus
    })
    
    wfun.pre <- function(predatmplus) {
      
      wpre.RI <- mplusObject(rdata = predatmplus, autov = T, 
                             DATA = "VARIANCES=NOCHECK;",
                             VARIABLE = "
 USEVARIABLES = trt pre scaled_w; 
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
 pre; pre on trt;
 
 %BETWEEN%
 pre; [pre];
 ",
   OUTPUT = "NOCHISQUARE TECH3;"
      )
      
      mplusname <- tempfile("wpre_RI", tmpdir = "tmp")
      fit.wpre.RI <- mplusModeler(wpre.RI, 
                                  dataout = paste0(mplusname, ".dat"),
                                  modelout = paste0(mplusname, ".inp"),
                                  run = 1) 
      
      est.wpre.RI <- fit.wpre.RI$results$parameters$unstandardized[1, c(3, 4, 6)]
      colnames(est.wpre.RI) <- c("Est","SE","pval")
      
      # pre_SL
      prepsw.design <- svydesign(ids=~school, nest=TRUE, 
                                 weights =~scaled_w, 
                                 data=predatmplus)
      est.wpre.SL <-summary(svyglm(formula = pre ~ trt, family=gaussian(),
                                   design=prepsw.design))$coefficients[2, c(1,2,4)]
      names(est.wpre.SL) <- c("Est","SE","pval")
      
      # combine model results
      data.frame(premod = c("preRI","preSL"), rbind(est.wpre.RI, est.wpre.SL))
      
    }
    
    
    wfun.pre(predatmplus.ipw) -> preest.ipw
    wfun.pre(predatmplus.clw) -> preest.clw
    
    ipw.preest <- data.frame(psmod=psname, weighting="IPW",
                             preest.ipw)
    
    clw.preest <- data.frame(psmod=psname, weighting="clustered",
                             preest.clw)
    
    rbind(clw.preest, ipw.preest)
    
  }
  lapply(1:length(psmod.List), estPre.fun) -> estpre.out.list
  
  
  # y model 
  
  wps <- 1
  estY.fun <- function(wps){
    psname <- names(wList)[[wps]]
    
    wList[[wps]] -> wdat
    
    datmplus.ipw <- with(wdat$wdimp, {
      datmplus <-data.frame(school, trt, y, 
                            scaled_w=scaled_ipw)
      datmplus
    })
    
    datmplus.clw <- with(wdat$wdatcl, {
      datmplus <-data.frame(school, trt, y, 
                            scaled_w=scaled_clw)
      datmplus
    })
    
    wfun.Y <- function(datmplus) {
      
      wY.RI <- mplusObject(rdata = datmplus, autov = T, 
                           DATA = "VARIANCES=NOCHECK;",
                           VARIABLE = "
                           USEVARIABLES = trt y scaled_w; 
                           MISSING ARE ALL (999);
                           
                           CLUSTER = school;
                           
                           WITHIN = trt; !not model trt to be latent cluster mean centered; 
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
      
      mplusname <- tempfile("wY_RI", tmpdir = "tmp")
      fit.wY.RI <- mplusModeler(wY.RI, dataout = paste0(mplusname, ".dat"),
                                modelout = paste0(mplusname, ".inp"), run = 1) # If "run" greater than one, the model is bootstrapped with run replications as well as the basic model.
      #fit.wY.RI <- readModels( paste0(  "wY_RI.out"), what = c("all"))
      
      est.wY.RI <- fit.wY.RI$results$parameters$unstandardized[1, c(3, 4, 6)]
      colnames(est.wY.RI) <- c("Est","SE","pval")
      
      # Y_SL
      psw.design <- svydesign(ids=~school, nest=TRUE, 
                              weights =~scaled_w, data=datmplus)
      est.wY.SL <-summary(svyglm(formula = y~trt, family=gaussian(), design=psw.design))$coefficients[2, c(1,2,4)]
      names(est.wY.SL) <- c("Est","SE","pval")
      
      data.frame(ymod = c("yRI","ySL"), rbind(est.wY.RI, est.wY.SL))
    }
    
    
    wfun.Y(datmplus.ipw) -> est.ipw
    wfun.Y(datmplus.clw) -> est.clw
    
    ipw.est <- data.frame(psmod=psname, weighting="IPW",
                          est.ipw)
    
    clw.est <- data.frame(psmod=psname, weighting="clustered",
                          est.clw)
    
    rbind(clw.est, ipw.est)
    
  }
  
  lapply(1:length(psmod.List), estY.fun) -> estY.out.list
  
  
  cbl.out.df <- do.call(rbind, cbl.out.list) 
  estpre.out.df <- do.call(rbind, estpre.out.list) 
  estY.out.df <- do.call(rbind, estY.out.list) 
  
  one.cbl.out <- data.frame(MI=MIname, cbl.out.df) 
  one.estpre.out <- data.frame(MI=MIname, estpre.out.df) 
  one.estY.out <- data.frame(MI=MIname, estY.out.df) 
  
  list(one.cbl.out=one.cbl.out,
       one.estpre.out=one.estpre.out, 
       one.estY.out=one.estY.out,
       psconv.List=psconv.List
  )
}

One.Imp(p=1, imputedlist = imputedlist.micex, MIname = "mice_x") -> micex1

View(micex1$one.cbl.out)


# parallel
micex.res.list <- mclapply(
  1:length(imputedlist.micex), FUN = One.Imp,
  imputedlist = imputedlist.micex, MIname = "mice_x",
  mc.preschedule = T, mc.cores = 20)

save(micex.res.list, file = "micex.res.list.RData")

load("micex.res.list.RData")


# Pool the results "micex" ----
p<-1

# covariate balance
pool.cbl.out <- lapply(1:length(micex.res.list), FUN = function(p){
  data.frame(imp=p, micex.res.list[[p]]$one.cbl.out)
}
)
pool.cbl <- do.call(rbind, pool.cbl.out)

micex.summ.cbl <- pool.cbl %>% 
  group_by(MI, psmod, weighting
           # , covname, imp
  ) %>%
  summarise_at(vars(w.smd, u.smd),
               list( max= ~max(abs(.)), 
                     # n05=~sum(abs(.)>.05), 
                     mean= ~mean(abs(.)) )
  ) %>% relocate(contains("u.smd"), .after = contains("w.smd"))


# BY level-2 covariates
micex.summ.cbl.L2L1 <- pool.cbl %>% 
  mutate( covLevel = ifelse(covname %in% c("bysctrl",  "byurban",  "byregion"), "L2", "L1")  ) %>%
  group_by(MI, psmod, weighting
           ,covLevel # , covname
  ) %>%
  summarise_at(vars(w.smd, u.smd),
               list( max= ~max(abs(.)), 
                    # n05=~sum(abs(.)>.05), 
                     mean= ~mean(abs(.)) )
  ) %>% relocate(contains("u.smd"), .after = contains("w.smd"))




# # estimates pre
# pool.estpre.out <- lapply(1:length(micex.res.list), FUN = function(p){
#   data.frame(imp=p, micex.res.list[[p]]$one.estpre.out)
# }
# )
# pool.estpre <- do.call(rbind, pool.estpre.out)
# 
# micex.estpre <- pool.estpre %>%
#   group_by(MI, psmod, weighting, premod) %>%
#   summarize(
#     poolval = testEstimates(qhat = Est, uhat = SE^2 )$estimates[1,c(1,2,5)],
#     name = c("Est","SE","pval")
#   ) %>%
#   pivot_wider(names_from = c("name"), values_from = c("poolval"))



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

predMat[ , !colnames(predMat)%in%c("school") ] <- 3 # 3 = overall + group-level effect
#predMat[ c("bygpared_College","bygrdrpt_yes" ) , !colnames(predMat)%in%c("school") ] <- 1
# 4 = individual-level (random) and group-level (fixed) effect
predMat[ ,'school'] <- -2 # -2 = cluster variable

diag(predMat) <-0

meth <- make.method(exdatmi)
meth[ c("bymathse","byenglse", "byactctl", "byconexp", "byinstmo","bystprep",
        "BYSES1" ) ] <- "2l.pan" 
#meth[ c("bygpared_College","bygrdrpt_yes") ] <- "2l.binary" 
#meth[ !colnames(predMat)%in%c("trt","y","school")  ] <- "2l.pmm" 
# "pmm works best with large samples, and provides imputations that possess many characteristics of the complete data. " (https://stefvanbuuren.name/fimd/sec-pmm.html); "for multi-category variables, but 2l.pmm may be a workable alternative." (Enders 2018 chapter); "PMM has several advantages..." (Leite et al., 2021)

# mice_RIx <- mice::futuremice(data = exdatmi, m = 4, n.core=4,
#            method = meth, predictorMatrix = predMat,
#            maxit = 2, parallelseed = 12345 )

mice_RIx <- mice( data = exdatmi, m = 20,
                  method = meth, predictorMatrix = predMat,
                  maxit = 50, seed = 12345 )
save(mice_RIx, file = "Eg_mice_RIx_50.RData")

#plot(mice_RIx)


load("Eg_mice_RIx_50.RData")

complete(mice_RIx, action = "all", include = FALSE) -> imputedlist.miceRIx

# compare miceRIx vs micex
mice_RIx$method
lmer(byactctl~1+(1|school), imputedlist.miceRIx$`1`)
lmer(byactctl~1+(1|school), imputedlist.micex$`1`) # smaller ICC


# Each imputed dataset MI_RE (same as One.Imp) ---------------------
p<-1
imputedlist <- imputedlist.miceRIx
MIname <- "mice_RIx"

#One.Imp(p=1, imputedlist = imputedlist.miceRIx, MIname = "mice_RIx") -> miceRIx1

#View(micex1$one.cbl.out)


# parallel
miceRIx.res.list <- mclapply(
  1:length(imputedlist.miceRIx), FUN = One.Imp,
  imputedlist = imputedlist.miceRIx, MIname = "mice_RIx",
  mc.preschedule = T, mc.cores = 20)

save(miceRIx.res.list, file = "miceRIx.res.list.RData")


load("~/Library/CloudStorage/Box-Box/Labs/PS_Missingdata/Example/miceRIx.res.list.RData")



# Pool the results "miceRIx" -----
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



# BY level-2 covariates
miceRIx.summ.cbl.L2L1 <- pool.cbl %>% 
  mutate( covLevel = ifelse(covname %in% c("bysctrl",  "byurban",  "byregion"), "L2", "L1")  ) %>%
  group_by(MI, psmod, weighting
           ,covLevel # , covname
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

# # estimates pre
# pool.estpre.out <- lapply(1:length(miceRIx.res.list), FUN = function(p){
#   data.frame(imp=p, miceRIx.res.list[[p]]$one.estpre.out)
# }
# )
# pool.estpre <- do.call(rbind, pool.estpre.out)
# 
# miceRIx.estpre <- pool.estpre %>%
#   group_by(MI, psmod, weighting, premod) %>%
#   summarize(
#     poolval = testEstimates(qhat = Est, uhat = SE^2 )$estimates[1,c(1,2,5)],
#     name = c("Est","SE","pval")
#             ) %>%
#   pivot_wider(names_from = c("name"), values_from = c("poolval"))





# Collect the results -----------------

# covariate balance
summ.cbl <- rbind(ccx.summ.cbl, micex.summ.cbl, miceRIx.summ.cbl)

summ.cbl.L2L1 <- rbind(ccx.summ.cbl.L2L1, micex.summ.cbl.L2L1, miceRIx.summ.cbl.L2L1)

estY <- rbind(ccx.estY, micex.estY, miceRIx.estY)
estY$Est_Y0minusY1 <- (-estY$Est)*10 # re-scale to the origianl standard deviation of 10
estY$SE_rescale <- estY$SE*10

write.csv(summ.cbl, "covariate balance CC MI.csv", row.names = F)
write.csv(summ.cbl.L2L1, "covariate balance-byLevel CC MI.csv", row.names = F)
write.csv(estY, "estimates pooled ATE.csv", row.names = F)







# Make sense of the results
mice_RIx$method
imputedlist.miceRIx$`1` -> miceRIx1
miceRIx1 %>% group_by(school) %>%
  summarise_all( ~mean(., na.rm=T) )->cmeans

cor(cmeans, use = "complete.obs")->clcor



