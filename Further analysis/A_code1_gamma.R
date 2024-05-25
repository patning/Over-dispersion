library("pomp")

run_level <- 2
measles_Np <-          1000
measles_Nmif <-        50
measles_Nreps_eval <-  10
measles_Nreps_local <- 20
measles_Nreps_global <-20
measles_Nsim <-        100
measles_profile_pts <-  5
measles_profile_Nreps <- 1

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
  filter(town=="London" & year>=1950 & year<1956) %>%
  mutate(
    time=(julian(date,origin=as.Date("1950-01-01")))/365.25+1950
  ) %>%
  filter(time>1950 & time<1956) %>%
  select(time,cases) -> dat

dat$cases=dat$cases+1
############################################################################
############################Model setup#####################################
############################################################################

rproc_gamma <- Csnippet("
  double beta, br, seas, foi, dw, births;
  double rate[6], trans[6];
  
  // cohort effect
  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt) 
    br = cohort*birthrate/dt + (1-cohort)*birthrate;
  else 
    br = (1.0-cohort)*birthrate;

  // term-time seasonality
  t = (t-floor(t))*365.25;
  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
    else
      seas = 1.0-amplitude;

  // transmission rate
  beta = R0*(gamma+mu)*seas;
  // expected force of infection
  foi = beta*pow(I+iota,alpha)/pop;
  // white noise (extrademographic stochasticity)
  dw = rgammawn(sigmaSE,dt);

  rate[0] = foi*dw/dt;  // stochastic force of infection
  rate[1] = mu;             // natural S death
  rate[2] = sigma;        // rate of ending of latent stage
  rate[3] = mu;             // natural E death
  rate[4] = gamma;        // recovery
  rate[5] = mu;             // natural I death

  // Poisson births
  births = rpois(br*dt);
  
  // transitions between classes
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,E,&rate[2],dt,&trans[2]);
  reulermultinom(2,I,&rate[4],dt,&trans[4]);

  S += births   - trans[0] - trans[1];
  E += trans[0] - trans[2] - trans[3];
  I += trans[2] - trans[4] - trans[5];
  R = pop - S - E - I;
  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
  C += trans[4];           // true incidence
")

## ----rinit-------------------------------------------------

rinit_gamma <- Csnippet("
  double m = pop/(S_0+E_0+I_0+R_0);
  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);
  W = 0;
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

#double tol = 1.0e-18;

dat %>%
  pomp(t0=with(dat,2*time[1]-time[2]), 
       times="time", 
       rprocess=euler(rproc_gamma,delta.t=1/365.25),  #euler.sim(rproc,delta.t=1/365.25)
       rinit=rinit_gamma, #initializer=initlz,
       dmeasure=dmeas,
       rmeasure=rmeas,
       covar=covariate_table(covar,times="time"),
       accumvars=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi","cohort","amplitude",
                    "S_0","E_0","I_0","R_0")
  ) -> measles_pomp_gamma

measles_params_guess_gamma<-c(R0=56.8, mu=0.02, sigma=28.9, gamma=30.4,alpha=0.976, iota=2.9, 
                             rho=0.488, psi=0.116, cohort=0.56, amplitude=0.55, 
                             S_0=0.0297, E_0=5.17e-05, I_0=5.14e-05, R_0=0.97, 
                             sigmaSE=0.088)

pt_gamma <- parameter_trans(
  log=c("sigma","gamma","sigmaSE","psi","R0","alpha", "iota"),
  logit=c("amplitude","cohort","rho"),
  barycentric=c("S_0","E_0","I_0","R_0")
)

dat %>% 
  pomp(t0=with(dat,2*time[1]-time[2]),
       time="time",
       params=measles_params_guess_gamma,
       rprocess=euler(rproc_gamma,delta.t=1/365.25),
       rinit=rinit_gamma,
       dmeasure=dmeas,
       rmeasure=rmeas,
       partrans=pt_gamma,
       covar=covariate_table(covar,times="time"),
       accumvars=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi","cohort","amplitude",
                    "S_0","E_0","I_0","R_0")
  ) -> measles_pomp_gamma
############################################################################
############The following is global MLE#####################################
############################################################################

measles_rw.sd_gamma = rw_sd(
  R0=0.02,sigma=0.02,gamma=0.02,
  psi=0.02,cohort=0.02,amplitude=0.02,
  sigmaSE=0.02, rho=0.02, alpha=0.005, iota=0.02,
  S_0=ivp(0.2),E_0=ivp(0.2),I_0=ivp(0.2),R_0=ivp(0.2))

estpars_gamma <- setdiff(names(measles_params_guess_gamma),
                         c("R0"))

measles_params_guess.t_gamma <- partrans(measles_pomp_gamma,measles_params_guess_gamma,"toEst")

measles_params_guess.t.hi_gamma <- measles_params_guess.t.lo_gamma <- measles_params_guess.t_gamma
measles_params_guess.t.lo_gamma[estpars_gamma] <- measles_params_guess.t_gamma[estpars_gamma]-log(2)
measles_params_guess.t.hi_gamma[estpars_gamma] <- measles_params_guess.t_gamma[estpars_gamma]+log(2)


profile_design(
  R0=seq(from=log(10),to=log(60),length=measles_profile_pts),
  lower=measles_params_guess.t.lo_gamma,
  upper=measles_params_guess.t.hi_gamma,
  nprof=measles_profile_Nreps
) -> pd_gamma

pd_gamma <- as.data.frame(t(partrans(measles_pomp_gamma,t(pd_gamma),"fromEst")))


library(doParallel)
registerDoParallel()
library(doRNG)
registerDoRNG(3899882)

bake(file=paste0(rds_files_dir,"A_global_gamma1.rds"),{
  
  foreach (
    p=iter(pd_gamma,"row"),
    .combine=bind_rows, .errorhandling="remove", .inorder=FALSE
  ) %dopar% {
    
    
    measles_pomp_gamma %>%
      mif2(
        params=p,
        Nmif = measles_Nmif,
        rw.sd = measles_rw.sd_gamma,
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
  filter(is.finite(loglik)) -> global_gamma

global_gamma[which.max(global_gamma[,17]),]

