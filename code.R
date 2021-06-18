library("pomp")

run_level <- 3
measles_Np <-          switch(run_level,100, 1e3, 5e3)
measles_Nmif <-        switch(run_level, 10, 100, 200)
measles_Nreps_eval <-  switch(run_level,  2,  10,  20)
measles_Nreps_local <- switch(run_level, 10,  20,  40)
measles_Nreps_global <-switch(run_level, 10,  20, 100)
measles_Nsim <-        switch(run_level, 50, 100, 500) 
measles_profile_pts <-  switch(run_level,  3,  5,  30)
measles_profile_Nreps <- switch(run_level, 2,  3,  20)

rds_files_dir <- paste0("rds_",run_level,"/")
if(!dir.exists(rds_files_dir)) dir.create(rds_files_dir)

library(readr)
library(tibble)
library(magrittr)
library(dplyr)


############################################################################
############################Load data#######################################
############################################################################

daturl <- "https://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path(tempdir(),"twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

demog %>%
  filter(town=="London") %>%
  select(-town) -> demogLondon

demogLondon %>%
  summarize(
    time=seq(from=min(year),to=max(year),by=1/12),
    pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
    birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
  ) -> covar

measles %>%
  mutate(year=as.integer(format(date,"%Y"))) %>% 
  filter(town=="London" & year>=1950 & year<1964) %>%
  mutate(
    time=(julian(date,origin=as.Date("1950-01-01")))/365.25+1950
  ) %>%
  filter(time>1950 & time<1964) %>%
  select(time,cases) -> dat

############################################################################
############################Model setup#####################################
############################################################################
rproc <- Csnippet("
  double beta, br, seas, foi,  births;
  double rate[6], trans[6];
  
  // cohort effect
  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
    br = cohort*birthrate/dt + (1-cohort)*birthrate;
  else
    br = (1.0-cohort)*birthrate;

  // term-time seasonality
  t = (t-floor(t))*365.25;
  if ((t>=7 && t<=100) ||
      (t>=115 && t<=199) ||
      (t>=252 && t<=300) ||
      (t>=308 && t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
  else
      seas = 1.0-amplitude;

  // transmission rate
  beta = R0*(gamma+mu)*seas;

  // expected force of infection
  foi = beta*pow(I+iota,alpha)/pop; 
  
  rate[0] = foi;        // rate of infection
  rate[1] = mu;         // natural S death
  rate[2] = sigma;      // rate of ending of latent stage
  rate[3] = mu;         // natural E death
  rate[4] = gamma;      // rate of recovery
  rate[5] = mu;         // natural I death

  // Poisson births
  births = rpois(br*dt);
  
  //the following block is newly added//
  double p_rproc, prob_SE_alpha, prob_SD_alpha, probs[3], Y_SE_gamma, Y_SD_gamma, Y_SS_gamma, Y_SUM_gamma;
  int vals[3];
  
  p_rproc = C_beta*(1-exp(-(rate[0]+rate[1])*dt))/(rate[0]+rate[1]);
  prob_SE_alpha = p_rproc*rate[0];
  prob_SD_alpha = p_rproc*rate[1];
  
  Y_SE_gamma=rgamma(prob_SE_alpha,1);
  Y_SD_gamma=rgamma(prob_SD_alpha,1);
  Y_SS_gamma=rgamma(C_beta-prob_SE_alpha-prob_SD_alpha,1);
  Y_SUM_gamma=Y_SE_gamma+Y_SD_gamma+Y_SS_gamma;

  probs[0]=Y_SE_gamma/Y_SUM_gamma;
  probs[1]=Y_SD_gamma/Y_SUM_gamma;
  probs[2]=1-probs[0]-probs[1];
  rmultinom(S, probs, 3, vals);
  trans[0] = vals[0];
  trans[1] = vals[1];
  //the above block is newly added//

  // transitions between classes
  reulermultinom(2,E,&rate[2],dt,&trans[2]);
  reulermultinom(2,I,&rate[4],dt,&trans[4]);

  S += births   - trans[0] - trans[1];
  E += trans[0] - trans[2] - trans[3];
  I += trans[2] - trans[4] - trans[5];
  R = pop - S - E - I;
  C += trans[4];           // true incidence
")

## ----rinit-------------------------------------------------
rinit <- Csnippet("
  double m = pop/(S_0+E_0+I_0+R_0);
  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);
  C = 0;
")

## ----dmeasure-------------------------------------------------
dmeas <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 0.0; 
  if (cases > 0.0) {
    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)
           - pnorm(cases-0.5,m,sqrt(v)+tol,1,0) + tol;
  } else {
    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0) + tol;
  }
  if (give_log) lik = log(lik);
")

## ----rmeasure-------------------------------------------------
rmeas <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 0.0;
  cases = rnorm(m,sqrt(v)+tol);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")


dat %>%
  pomp(t0=with(dat,2*time[1]-time[2]), 
       times="time", 
       rprocess=euler(rproc,delta.t=1/365.25),  
       rinit=rinit,
       dmeasure=dmeas,
       rmeasure=rmeas,
       covar=covariate_table(covar,times="time"),
       accumvars=c("C"),
       statenames=c("S","E","I","R","C"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","C_beta", "psi","cohort","amplitude",
                    "S_0","E_0","I_0","R_0") 
  ) -> measles_pomp


measles_params_guess<-c(R0=46.80919, mu=0.02, sigma=24.30255, gamma=46.26716,alpha=0.976, iota=2.9, 
                        rho=0.488, psi=0.1198564, cohort=0.988439, amplitude=0.5043777, 
                        S_0=0.03184962, E_0=8.416195e-05, I_0=2.449359e-05, R_0=0.9680417, 
                        C_beta=739.9351)


pt <- parameter_trans(
  log=c("sigma","gamma","C_beta","psi","R0","alpha", "iota"),
  logit=c("cohort","amplitude","rho"),
  barycentric=c("S_0","E_0","I_0","R_0")
) 

dat %>% 
  pomp(t0=with(dat,2*time[1]-time[2]),
       time="time",
       params=measles_params_guess,
       rprocess=euler(rproc,delta.t=1/365.25),
       rinit=rinit,
       dmeasure=dmeas,
       rmeasure=rmeas,
       partrans=pt,
       covar=covariate_table(covar,times="time"),
       accumvars=c("C"),
       statenames=c("S","E","I","R","C"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","C_beta","psi","cohort","amplitude",
                    "S_0","E_0","I_0","R_0")
  ) -> measles_pomp



  
############################################################################
############The following is global MLE#####################################
############################################################################
measles_rw.sd = rw.sd(
  R0=0.02,sigma=0.02,gamma=0.02,
  psi=0.02,cohort=0.02,amplitude=0.02,
  C_beta=0.02, rho=0.02, alpha=0.005, iota=0.02,
  S_0=ivp(0.2),E_0=ivp(0.2),I_0=ivp(0.2),R_0=ivp(0.2))

estpars <- setdiff(names(measles_params_guess),
                   c("R0","mu"))

measles_params_guess.t <- partrans(measles_pomp,measles_params_guess,"toEst")

measles_params_guess.t.hi <- measles_params_guess.t.lo <- measles_params_guess.t
measles_params_guess.t.lo[estpars] <- measles_params_guess.t[estpars]-log(2)
measles_params_guess.t.hi[estpars] <- measles_params_guess.t[estpars]+log(2)


profile_design(
  R0=seq(from=log(10),to=log(60),length=measles_profile_pts),
  lower=measles_params_guess.t.lo,
  upper=measles_params_guess.t.hi,
  nprof=measles_profile_Nreps
) -> pd

pd <- as.data.frame(t(partrans(measles_pomp,t(pd),"fromEst")))


library(doParallel)
registerDoParallel()
library(doRNG)
registerDoRNG(3899882)

bake(file=paste0(rds_files_dir,"global-profile.rds"),{
  
  foreach (
    p=iter(pd,"row"),
    .combine=bind_rows, .errorhandling="remove", .inorder=FALSE
  ) %dopar% {
    
    
    measles_pomp %>%
      mif2(
        params=p,
        Nmif = measles_Nmif,
        rw.sd = measles_rw.sd,
        Np = measles_Np,
        cooling.type = "geometric",
        cooling.fraction.50 = 0.1
      ) %>%
      mif2() -> mf
    
    pf <- replicate(measles_Nreps_eval, pfilter(mf, Np = measles_Np))
    ll <- sapply(pf,logLik)
    ll <- logmeanexp(ll, se = TRUE)
    
    data.frame(
      as.list(coef(mf)),
      loglik = ll[1],
      loglik.se = ll[2]
    )
  }
}) %>%
  filter(is.finite(loglik)) -> global_prof


global_prof[which.max(global_prof[,17]),]
#               R0     R0.1   mu    sigma    gamma    alpha     iota       rho
# se...316 34.08814 34.08814 0.02 52.71022 22.88241 1.017399 55.08132 0.4915295
#             psi    cohort    amplitude        S_0          E_0          I_0
# se...316 0.1184038 0.9998928 0.4802633 0.03197983 6.989526e-05 4.521218e-05
#               R_0  C_beta    loglik loglik.se
# se...316 0.9679051 652.763 -3803.166 0.4209594

############################################################################
############The following is R0 profile#####################################
############################################################################
measles_rw.sd = rw.sd(
  sigma=0.02,gamma=0.02,
  psi=0.02,cohort=0.02,amplitude=0.02,
  C_beta=0.02, rho=0.02, iota=0.02,
  S_0=ivp(0.2),E_0=ivp(0.2),I_0=ivp(0.2),R_0=ivp(0.2))

estpars <- setdiff(names(measles_params_guess),
                   c("R0","mu","alpha"))

measles_params_guess["alpha"] <- 1
measles_params_guess.t <- partrans(measles_pomp,measles_params_guess,"toEst")

measles_params_guess.t.hi <- measles_params_guess.t.lo <- measles_params_guess.t
measles_params_guess.t.lo[estpars] <- measles_params_guess.t[estpars]-log(2)
measles_params_guess.t.hi[estpars] <- measles_params_guess.t[estpars]+log(2)


profile_design(
  R0=seq(from=log(10),to=log(60),length=measles_profile_pts),
  lower=measles_params_guess.t.lo,
  upper=measles_params_guess.t.hi,
  nprof=measles_profile_Nreps
) -> pd

pd <- as.data.frame(t(partrans(measles_pomp,t(pd),"fromEst")))



library(doParallel)
registerDoParallel()
library(doRNG)
registerDoRNG(3899882)

bake(file=paste0(rds_files_dir,"R0-profile.rds"),{
  
  foreach (
    p=iter(pd,"row"),
    .combine=bind_rows, .errorhandling="remove", .inorder=FALSE
  ) %dopar% {
    
    
    measles_pomp %>%
      mif2(
        params=p,
        Nmif = measles_Nmif,
        rw.sd = measles_rw.sd,
        Np = measles_Np,
        cooling.type = "geometric",
        cooling.fraction.50 = 0.1
      ) %>%
      mif2() -> mf
    
    pf <- replicate(measles_Nreps_eval, pfilter(mf, Np = measles_Np))
    ll <- sapply(pf,logLik)
    ll <- logmeanexp(ll, se = TRUE)
    
    data.frame(
      as.list(coef(mf)),
      loglik = ll[1],
      loglik.se = ll[2]
    )
  }
}) %>%
  filter(is.finite(loglik)) -> R0_prof


############################################################################
############The following is R0 profile plot################################
############################################################################

R0_prof %>%
  mutate(R0=exp(signif(log(R0),5))) %>%
  group_by(R0) %>%
  filter(rank(-loglik)<=1) %>%
  ungroup() -> R0_prof

line.color <- 'gray20'
x<-R0_prof

par(mfrow=c(1,1), cex.lab=1.5, cex.axis=1)
plot.new()
mtext(side=2,line=2.8,text="profile log likelihood",adj=0.5)
mtext(side=1,line=2.5,text=expression(R[0]),adj=0.5)

fit <- loess(loglik~R0,data=x,span=0.8)
nd <- data.frame(R0=with(x,seq(from=min(R0),to=max(R0),length=100)))
nd$loglik <- predict(fit,newdata=nd)
cutoff <- max(nd$loglik,na.rm=T)-0.5*qchisq(p=0.95,df=1)
ci.lond <- conf.int <- range(nd$R0[nd$loglik>cutoff],na.rm=T)
par(fig=c(0,1,0,1),mar=c(4,4,2,2),new=T)
plot(
  loglik~R0,
  data=x,
  font=2,
  bty='l',
  ann=F,
  xaxt='n',
  xlim=c(20,70),
  ylim=max(x$loglik)+c(-5,1)
)
axis(side=1,at=seq(20,70,by=10))
lines(loglik~R0,data=nd,col=line.color)
abline(h=cutoff,lty='33')
abline(v=conf.int,lty='63')
conf.int
#[1] 31.21107 47.37178
