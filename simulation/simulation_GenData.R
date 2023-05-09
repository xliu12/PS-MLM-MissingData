
library(mvtnorm)

GenData <- function( seed=9999,
                     J=100, # number of clusters 
                     nj=25, # cluster size 
                     x_missing_rate=0.20, # missing rate of each individual covariate 
                     # unobserved L2 moderator/confounder influence missing status
                     heter_treat =  "treat_by_z", # "treat_homo"
                     x_on_z =0, #0.5 # effect of cluster-level confounder Z on individual covariate
                     missingness= "MCAR", # "MAR_z"
                     iccx = 0.1, # 0.5
                     b_xcontextual = 0, #1 contextual effect of X1 on Y, i.e., effect of X1b, conditional between-cluster component of X1 conditional on Z, on outcome y
                     a_xcontextual = 0, #1 contextual effect of X1 on latent treatment `trt_ast`
                     ATE = 0.5 # true treatment effect
) {
  # seed=111
  set.seed(seed = seed)
  
  # J <- 100   
  # nj <- 100  
  
  N <- J*nj
  school <- rep(1:J, each=nj)
  
  #### covariates ----
  
  num_x <- 5;
  num_xmiss <- 4  # incomplete X
  num_z <- 1 # unobserved Z
  
  z <- rnorm(J, 0, 1)[ rep(1:J, each=nj)]
  
  # correlation between incomplete covariates' within-cluster errors 
  corr_ex <- 0.1  # covariates should cover heterogeneous domains
  Cov_ex <- diag( 1, nrow = num_x ); Cov_ex[ lower.tri(Cov_ex)]= Cov_ex[upper.tri(Cov_ex)] <-(corr_ex)
  ex <- rmvnorm(N, sigma = Cov_ex )
  
  x5 <- ex[,5]
     
  xb_ast  <- rmvnorm(J, sigma = diag(1*iccx/(1-iccx), nrow = num_xmiss) )[ rep(1:J, each=nj), ]
  
  x_ast <- xb_ast[,1:4] + x_on_z* z + ex[,1:4]
  
  x <- cbind(x_ast, x5) # cor(x)
  colnames(x)<-paste0( 'x.', 1:num_x)
  
  # missing indicators -----
  
  # for each of x1,x2,x3,x4 that have missing values
  num_xmiss <- 4; 
  
  if( missingness == "MCAR" ){
    
    g_zobs<-0    
  }
  if( missingness == "MAR_z"){
    
    g_zobs <- (-1)*sqrt(0.35) 
  }
  
  g0 <- qnorm(x_missing_rate,mean=0, sd = sqrt(g_zobs^2+1), lower.tail = TRUE); 
  rx <- g0 + g_zobs*z + replicate(num_xmiss, rnorm(N))  
  
  Rx <- ifelse(rx>0, 1, 0)
  
  x_obs <- x
  x_obs[,1:num_xmiss] <- ifelse( rx>0, NA, x[, 1:num_xmiss] )
  
  
  #--- treatment ----
  icct <- .2  
  a_0<-0
  a_x1 <-c( sqrt( ((0.26)/(1-0.26))*(1-icct)/sum(Cov_ex) ))
  a_z1 <- c(sqrt( ((0.26)/(1-0.26))*(icct)/1 ))
  if(iccx !=0){
    a_xc1 <- sqrt( 2* ((0.26)/(1-0.26))*(icct) / (1*(iccx)/(1-iccx))   ) 
  }
  if(iccx ==0){
    a_xc1 <- 0
  }
  et <- rlogis(N) * sqrt(1-icct) 
  sd_logit <- sqrt(pi^2/3)
  ut <- rnorm( J, sd = sqrt(icct) )[rep(1:J,each=nj) ]
  
  t_ast <- a_0+ sd_logit* (rowSums(a_x1*x) + a_z1*z) + 
    sd_logit*ut + et  + sd_logit*(a_xcontextual* a_xc1*xb_ast[,1] )  
  
  trt <- as.numeric(t_ast > 0)
  
  #### Outcome ----
  
  iccy <- .2  # ICC of Y0.residual
  b_0 <- 0
  b_x1 <-sqrt( ((0.26)/(1-0.26))*(1-iccy)/sum(Cov_ex) ) 
  b_z1 <- sqrt( ((0.26)/(1-0.26))*(iccy)/1 )
  
  # ATE <- 0.5 
  tau0 <- ATE
  
  if( heter_treat == "treat_homo" ){
    tau_z1 <- 0
    tau_x1 <-0
  }
  
  if( heter_treat == "treat_by_z" ){
    
    tau_z1 <- sqrt( 0.5*(iccy)/1 )
    tau_x1 <-0 
  }
  
  #---- outcome
  uy <- rnorm( J, sd = sqrt(iccy) )[rep(1:J,each=nj)]
  ey <-  rnorm( N, sd=sqrt(1-iccy) )
  vy <- rnorm( J, sd = sqrt( iccy) )[rep(1:J,each=nj)]
  
  # context effects of incomplete x[,1] e.g., pretest
  if(iccx!=0){
    b_xc1 <- sqrt( 2* ((0.26)/(1-0.26))*(iccy) / (1*(iccx)/(1-iccx))   ) # xb variance: (iccx)/(1-iccx)
    
  }
  if(iccx==0){
    b_xc1 <- 0
  }
  
  
  y <- b_0+ (tau0 +  tau_x1* x[,5] + tau_z1*z + vy)*trt + 
    rowSums( b_x1*x[,1:4] ) + b_x1*x[,5] + b_z1*z +
    (b_xcontextual* b_xc1*xb_ast[,1] ) + uy  + ey
  
  
  data.frame(school= rep(1:J,each=nj), x=x, 
             trt=trt, y=y, x_obs=x_obs, Rx=Rx ) ->dat
  colnames(dat) <- c("school", paste0("x.", 1:num_x), "trt","y",paste0("x_obs.",1:num_x),
                      paste0("Rx.", 1:num_xmiss) )
  # lm(y~trt, dat)
  out <- mget( ls(), envir = environment() )
  
  return(out)
}

