library(mice)
library(miceadds)
library(mitml)
library(lme4)
library(MplusAutomation)
library(merDeriv)
library(sandwich)
library(survey)

source("simulation_GenData.R")
out_GenData <- GenData()

# naive: no covariate adjustment 
unadj.Outcome <- function(out.GenData = out_GenData ){
  datmplus <- out.GenData$dat[ ,c('school','trt','y'  )]
  
  
  Y_lm <- summary(lm( y~trt, data = datmplus ))
  
  out <- data.frame(Ymodel = "ySL", Estimate = NA, Std.Error=NA, pval=NA,
                    psconv=1, PSmodel="no")
  
  out[1, 2:4] <- as.numeric(Y_lm$coefficients[2, c(1,2,4)])
  
  return(out)
}


# PS, outcome model producing ATE
IPW.Outcome <- function( impslist = imps, # imps <- as.mitml.list(list( out.GenData$dat))
                         psmodel = "ps_REx_benchmark", # "benchmark": no-missingness
                         weighting="IPW", #"CLW", 
                         ymodel=c("ySL","yRE"), 
                         out.GenData = out_GenData
                         ) {
  
  clusname<-"school"
  zname<-"trt"; yname<-"y"
  
  num_x <- out.GenData$num_x
  
  ps_RExcmx <- paste0('trt ~ 1+ ', paste0( 'x_obs.', 1:num_x, collapse = '+'), '+',
                      paste0( 'cmx_obs.', 1:num_x, collapse = '+'),  '+ ( 1 | school )' )
  ps_RExcmx_benchmark <- paste0('trt ~ 1+ ', paste0( 'x.', 1:num_x, collapse = '+'), '+',
                                paste0( 'cmx.', 1:num_x, collapse = '+'),  '+ ( 1 | school )' )
  
  ps_REx <- paste0('trt ~ 1+ ', paste0( 'x_obs.', 1:num_x, collapse = '+'),  '+ ( 1 | school )' )
  ps_REx_benchmark <- paste0('trt ~ 1+ ', paste0( 'x.', 1:num_x, collapse = '+'),  '+ ( 1 | school )' )
  
  ps_FEx <- paste0('trt ~ 1+ ', paste0( 'x_obs.', 1:num_x, collapse = '+'), '+ as.factor(school) ' )
  ps_FEx_benchmark <- paste0('trt ~ 1+ ', paste0( 'x.', 1:num_x, collapse = '+'), '+ as.factor(school) ' )
  
  
  # with missing indicators
  ps_REx_mp <- paste0('trt ~ 1+ ', paste0( 'x_obs.', 1:num_x, collapse = '+'), '+',
                      paste0('Rx.', 1:out.GenData$num_xmiss, collapse = '+' ),
                      '+ ( 1 | school )' )
  ps_FEx_mp <- paste0('trt ~ 1+ ', paste0( 'x_obs.', 1:num_x, collapse = '+'), '+',
                      paste0( 'Rx.', 1:out.GenData$num_xmiss, collapse = '+'),
                      '+ as.factor(school) ')
  
  
  ipw_imps <- with(
    impslist, 
    {
      # ps methods ----
      
      out_s <-NULL
      for ( pm in psmodel ){ # pm<-psmodel[1]
        psmodel_one <- pm
        # with x, x-cluster-means, Model can fail to converge
        if( psmodel_one=="ps_RExcmx_benchmark" ){ 
          cmx <- aggregate( data.frame(x.1, x.2, x.3, x.4, x.5), by=list(school), FUN=mean, na.rm=TRUE  )
          colnames(cmx) <- c("school", paste0("cmx.", 1:5 ) )
          dxcmx <- merge(data.frame( school, x.1, x.2, x.3, x.4, x.5, trt, y ), cmx, by ="school" )
          ps_RExcmx_benchmark<-paste0('trt ~ 1+ ', paste0( 'x.', 1:num_x, collapse = '+'), '+',
                                      paste0( 'cmx.', 1:num_x, collapse = '+'),  '+ ( 1 | school )' )
          
          psmod <- glmer( formula(ps_RExcmx_benchmark), family="binomial", data = dxcmx  ) 
          psconv <- as.numeric(!any( grepl("failed to converge", psmod@optinfo$conv$lme4$messages) ))
        }
        
        
        if( psmodel_one=="ps_RExcmx" ){ 
          cmx <- aggregate( data.frame(x_obs.1, x_obs.2, x_obs.3, x_obs.4, x_obs.5), by=list(school), FUN=mean, na.rm=TRUE  )
          colnames(cmx) <- c("school", paste0("cmx_obs.", 1:5 ) )
          dxcmx <- merge(data.frame( school, x_obs.1, x_obs.2, x_obs.3, x_obs.4, x_obs.5, trt, y ), cmx, by ="school" )
          psmod <- glmer( ps_RExcmx, family="binomial", data = dxcmx  ) 
          psconv <- !any( grepl("failed to converge", psmod@optinfo$conv$lme4$messages) )
          
        }
        
        
        if( psmodel_one=="ps_REx" ){ 
          
          psmod <- glmer( ps_REx, family="binomial"  ) 
          psconv <- !any( grepl("failed to converge", psmod@optinfo$conv$lme4$messages) )
        }
        
        if( psmodel_one=="ps_REx_benchmark" ){ 
          psmod <- glmer( ps_REx_benchmark, family="binomial"  ) 
          psconv <- !any( grepl("failed to converge", psmod@optinfo$conv$lme4$messages) )
        }
        
        
        
        if( psmodel_one=="ps_FEx" ){ 
          psmod <- glm( formula(ps_FEx), family="binomial") 
          psconv <- (psmod$converged)
        }
        if( psmodel_one=="ps_FEx_benchmark" ){ 
          psmod <- glm( formula(ps_FEx_benchmark), family="binomial") 
          psconv <- (psmod$converged)
        }
        
        # missing pattern indicators (Rx.1,...,Rx.4)
        if( psmodel_one=="ps_REx_mp" ){ 
          datmp <- data.frame( school, x_obs.1, x_obs.2, x_obs.3, x_obs.4, x_obs.5, trt, y,
                               Rx=out.GenData$Rx )
          psmod <- glmer( ps_REx_mp, family="binomial", data = datmp  ) 
          psconv <- !any( grepl("failed to converge", psmod@optinfo$conv$lme4$messages) )
        }
        
        if( psmodel_one=="ps_FEx_mp" ){ 
          datmp <- data.frame( school, x_obs.1, x_obs.2, x_obs.3, x_obs.4, x_obs.5, trt, y,
                               Rx=out.GenData$Rx )
          psmod <- glm( formula(ps_FEx_mp ), family="binomial", data = datmp  ) 
          psconv <- (psmod$converged)
          
        }
        
        
        
        ps<-psest <- fitted(psmod)
        w <- trt * (1/ps) + (1-trt) * (1/(1-ps))
        
        if(weighting=="IPW"){
          scaled_w <- (trt==1)*w /sum( w[trt ==1]) +  (trt==0)*w /sum( w[trt ==0])
        }
        
        if( weighting=="CLW" ){
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
        # 
        if( weighting!="CLW"){
          datmplus <-data.frame( school, trt, y, scaled_w)
    
        }
        
        if( weighting=="CLW" ){
          datmplus <-data.frame( school=datcl$school, trt=datcl$trt, y=datcl$y, scaled_w=datcl$cl_w)
          
        }
        
        
        out <- data.frame(Ymodel = ymodel, Estimate = NA, Std.Error=NA, pval=NA,
                          psconv=psconv, PSmodel=psmodel_one)
        
        # mplus ----
        if( length(grep("yR", ymodel))!=0 ){
          pid <- Sys.getpid()
          write.table(datmplus, paste0(pid,w[1],"datmplus.txt"), 
                      na="999", quote = F,row.names = F,col.names = F )
        }
        
        if( "yRE"%in% ymodel ) {
          cat(
            sep = " ",
            " Data: 
            file =", paste0(pid,w[1],"datmplus.txt" ) ,";
            VARIANCES=NOCHECK;
            
            Variable:
            NAMES = ", colnames(datmplus), ";
            USEVARIABLES = trt y scaled_w; 
            MISSING ARE ALL (999);
            CLUSTER = school;
            
            WITHIN = trt; 
            BETWEEN = ;
            WEIGHT = scaled_w; WTSCALE = unscaled;
            
            Analysis:
            Type = Twolevel Random;
            Estimator = MLR;
            
            Model:
            %WITHIN%
            y; y on trt;
            
            %BETWEEN%
            y; [y];
            
            Output:
            NOCHISQUARE TECH1 TECH3;
            ", 
            file = paste0(pid,w[1], "wY_RI.inp")
          )
          runModels( paste0(pid,w[1], "wY_RI.inp"), logFile = NULL)
          wY_RI_out=readModels( paste0(pid,w[1], "wY_RI.out"), what = c("all"))
          est_RI <- wY_RI_out$parameters$unstandardized[1, ]
          
          an.error.occured <- FALSE
          tryCatch( { est_RI[ ,c(3, 4, 6)] }
                    , error = function(e) {an.error.occured <<- TRUE})
          if( an.error.occured==FALSE) {
            out[out$Ymodel=="yRE", 2:4] <- est_RI[ ,c(3, 4, 6)] }
        }
        
        
        # end mplus 
        
        
        # single-level ----
        if("ySL" %in% ymodel){
          
          fit.wy <- lm( y~trt, data = datmplus, weights = scaled_w )
          est_tau <- fit.wy$coefficients["trt"] 
          out[out$Ymodel=="ySL", 2:4] <-c(est_tau, NA, NA) 
          
        }
        
        out_s <- data.frame(rbind(out_s, out ))
      }
      
      out_s
      #out=mget( ls(), envir = environment())
    } )
  
  
  
  return(ipw_imps)
}

# No-Missing
benchmark <- function(
    out.GenData=out_GenData, psmodel = "ps_FEx_benchmark", weighting="CLW",
    ymodel=c("yRE")
){
  imps <- as.mitml.list(list( out.GenData$dat))
  out_GenData <- out.GenData
  
  df<-NULL
  if("CLW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel,ymodel = ymodel,
                            weighting = "CLW" , out.GenData = out_GenData )
    out <- do.call(rbind, ipw_imps )
    
    df_wcl <- out; 
    df_wcl$weighting <-"CLW"
    
    df<-rbind( df, df_wcl)
  }
  if("IPW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel,ymodel = ymodel,
                            weighting = "IPW" , out.GenData = out_GenData )
    out <- do.call(rbind, ipw_imps )
    
    df_ipw <- out; 
    df_ipw$weighting <-"IPW"
    
    df<-rbind( df, df_ipw)
  }
  
  
  return(df)
}


complete_case <- function( out.GenData=out_GenData, psmodel = "ps_FEx",weighting = "CLW", ymodel=c("yRE") ){
  num_x <- out.GenData$num_x
  datx <- out.GenData$dat[ ,c('school','trt','y',paste0('x_obs.',1:num_x) )]
  dat_completex <- datx[ complete.cases(datx), ]
  imps <- as.mitml.list(list( dat_completex ))
  out_GenData <- out.GenData
  
  df<-NULL
  if("CLW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel,ymodel = ymodel,
                            weighting = "CLW" , out.GenData = out_GenData )
    out <- do.call(rbind, ipw_imps )
    
    df_wcl <- out; 
    df_wcl$weighting <-"CLW"
    
    df<-rbind( df, df_wcl)
  }
  if("IPW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel,ymodel = ymodel,
                            weighting = "IPW" , out.GenData = out_GenData )
    out <- do.call(rbind, ipw_imps )
    
    df_ipw <- out; 
    df_ipw$weighting <-"IPW"
    
    df<-rbind( df, df_ipw)
  }
  
  
  return(df)
}


### imputating  ----
num_imps <-20
burn <- 500 
miceseed <- 83492

# ignore clusters in the imputation model 
MI_SL <- function( out.GenData=out_GenData, psmodel = c("ps_FEx"), 
                    weighting = c("CLW"), ymodel=c("yRE") ) {
  num_x <- out.GenData$num_x
  datmice <- out.GenData$dat[ ,c('school','trt','y',paste0('x_obs.',1:num_x) )]
  
  meth <- make.method(datmice)
  meth[paste0('x_obs.',1:4) ] <- "pmm"
  predMat <- make.predictorMatrix(data = datmice)
  predMat[ ,'school'] <- 0
  
  mice_x <- mice( data = datmice, m = num_imps, method = meth, predictorMatrix = predMat,
                  maxit = burn, seed = miceseed )
  complete(mice_x, action = "all", include = FALSE) -> imputedlist
  
  imps <- as.mitml.list( imputedlist )
  
  out_GenData <- out.GenData
  df<-NULL
  if("CLW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel ,ymodel = ymodel,
                            weighting = "CLW" ,  out.GenData = out_GenData )
    
    poolest <- ipw_imps[[1]]
    poolest[, c('Estimate','Std.Error')]<-NA
    for(i in 1:nrow(ipw_imps[[1]])){
      oneparam_qu <- do.call(rbind, lapply(ipw_imps, function(oneipw_imps){
        oneipw_imps[i, c('Estimate','Std.Error')]
      }) )
      poolerror<-FALSE
      tryCatch({
        qhats <- oneparam_qu[complete.cases(oneparam_qu), 'Estimate']
        sehats <- oneparam_qu[complete.cases(oneparam_qu), 'Std.Error']
        uhats <- sehats^2
        poolone <- testEstimates(qhat = qhats, uhat = uhats)
        poolest[i, c('Estimate','Std.Error','pval')] <- poolone$estimates[ ,c('Estimate','Std.Error','P(>|t|)')]
      }, error = function(e) {poolerror <<- TRUE} )
      #print(poolest[i, ])
    }
    
    
    df <- poolest; 
    df$weighting <-"CLW"
  }
  if("IPW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel ,ymodel = ymodel,
                            weighting = "IPW" ,  out.GenData = out_GenData )
    
    poolest <- ipw_imps[[1]]
    poolest[, c('Estimate','Std.Error')]<-NA
    for(i in 1:nrow(ipw_imps[[1]])){
      oneparam_qu <- do.call(rbind, lapply(ipw_imps, function(oneipw_imps){
        oneipw_imps[i, c('Estimate','Std.Error')]
      }) )
      poolerror<-FALSE
      tryCatch({
        qhats <- oneparam_qu[complete.cases(oneparam_qu), 'Estimate']
        sehats <- oneparam_qu[complete.cases(oneparam_qu), 'Std.Error']
        uhats <- sehats^2
        poolone <- testEstimates(qhat = qhats, uhat = uhats)
        poolest[i, c('Estimate','Std.Error','pval')] <- poolone$estimates[ ,c('Estimate','Std.Error','P(>|t|)')]
      }, error = function(e) {poolerror <<- TRUE} )
      #print(poolest[i, ])
    }
    
    df <- poolest; 
    df$weighting <-"IPW"
    
  }
  
  rownames(df) <- NULL
  
  return(df)
}


# cluster random intercepts in the imputation model 
MI_RE <- function( out.GenData=out_GenData, psmodel = "ps_FEx" ,weighting = "CLW" ,  ymodel="yRE") {
  num_x <- out.GenData$num_x
  datmice <- out.GenData$dat[ ,c('school','trt','y',paste0('x_obs.',1:num_x) )]
  
  predMat <- make.predictorMatrix(data=datmice) # same as mice_ini$predictorMatrix
  
  predMat[ , paste0('x_obs.',1:num_x) ] <- 3 # 3 = overall + group-level effect
  # 4 = individual-level (random) and group-level (fixed) effect
  predMat[ , c('y') ] <- 3 ; 
  predMat[ , c('trt') ] <- 3 
  predMat[ ,'school'] <- -2 # -2 = cluster variable
  diag(predMat) <-0
  
  meth <- make.method(datmice)
  meth[paste0('x_obs.',1:4) ] <- "2l.pan"
  
  mice_RIx <- mice( data = datmice, m = num_imps, method = meth, predictorMatrix = predMat,
                    maxit = burn, seed = miceseed ) 
  complete(mice_RIx, action = "all", include = FALSE) -> imputedlist
  
  imps <- as.mitml.list( imputedlist )
  
  out_GenData <- out.GenData
  df<-NULL
  if("CLW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel ,ymodel = ymodel,
                            weighting = "CLW" ,  out.GenData = out_GenData )
    
    poolest <- ipw_imps[[1]]
    poolest[, c('Estimate','Std.Error')]<-NA
    for(i in 1:nrow(ipw_imps[[1]])){
      oneparam_qu <- do.call(rbind, lapply(ipw_imps, function(oneipw_imps){
        oneipw_imps[i, c('Estimate','Std.Error')]
      }) )
      poolerror<-FALSE
      tryCatch({
        qhats <- oneparam_qu[complete.cases(oneparam_qu), 'Estimate']
        sehats <- oneparam_qu[complete.cases(oneparam_qu), 'Std.Error']
        uhats <- sehats^2
        poolone <- testEstimates(qhat = qhats, uhat = uhats)
        poolest[i, c('Estimate','Std.Error','pval')] <- poolone$estimates[ ,c('Estimate','Std.Error','P(>|t|)')]
      }, error = function(e) {poolerror <<- TRUE} )
      #print(poolest[i, ])
    }
    
    df <- poolest; 
    df$weighting <-"CLW"
    
  }
  if("IPW" %in% weighting){
    ipw_imps <- IPW.Outcome(impslist = imps, psmodel = psmodel ,ymodel = ymodel,
                            weighting = "IPW" ,  out.GenData = out_GenData )
    
    poolest <- ipw_imps[[1]]
    poolest[, c('Estimate','Std.Error')]<-NA
    for(i in 1:nrow(ipw_imps[[1]])){
      oneparam_qu <- do.call(rbind, lapply(ipw_imps, function(oneipw_imps){
        oneipw_imps[i, c('Estimate','Std.Error')]
      }) )
      poolerror<-FALSE
      tryCatch({
        qhats <- oneparam_qu[complete.cases(oneparam_qu), 'Estimate']
        sehats <- oneparam_qu[complete.cases(oneparam_qu), 'Std.Error']
        uhats <- sehats^2
        poolone <- testEstimates(qhat = qhats, uhat = uhats)
        poolest[i, c('Estimate','Std.Error','pval')] <- poolone$estimates[ ,c('Estimate','Std.Error','P(>|t|)')]
      }, error = function(e) {poolerror <<- TRUE} )
      #print(poolest[i, ])
    }
    
    df <- poolest; 
    df$weighting <-"IPW"
}
  
  rownames(df) <- NULL
  
  return(df)
}





