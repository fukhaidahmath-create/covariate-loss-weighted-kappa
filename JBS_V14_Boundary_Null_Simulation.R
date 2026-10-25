# =============================================================================
# JBS V14 - TARGETED BOUNDARY-DOMINANCE VALIDATION (standalone supplement)
# Research diagnostic; DOES NOT alter or overwrite the original V13.3/V13.4 run.
# Works with base R. For cross-fitted mode, requires package brglm2.
# Two exact boundary DGM truth cases with P{h(X)=0}=0.35:
#   lower_zero: lower B = 0, upper B > 0;
#   upper_zero: lower B < 0, upper B = 0 (test labels interchanged).
# These are new stress-test scenarios, not part of the original 120 cells.
# The default 'oracle' analysis tests inferential geometry at KNOWN p,q.
# Set JBS_V14_ESTIMATOR='crossfit' to run the original logistic nuisance recipe.
# No claim of uniform validity arises from success in these simulations.
# =============================================================================

v14_options <- list(
  out_dir = Sys.getenv('JBS_V14_DIR', unset = file.path(getwd(), 'JBS_V14_BOUNDARY_RESULTS')),
  estimator = Sys.getenv('JBS_V14_ESTIMATOR', unset='oracle'),
  reps = as.integer(Sys.getenv('JBS_V14_REPS', unset='200')),
  ns = as.integer(strsplit(Sys.getenv('JBS_V14_NS', unset='250,858,1500'), ',', fixed=TRUE)[[1L]]),
  rhos = as.numeric(strsplit(Sys.getenv('JBS_V14_RHOS', unset='0.25,0.55,0.85'), ',', fixed=TRUE)[[1L]]),
  taus = as.numeric(strsplit(Sys.getenv('JBS_V14_TAUS', unset='0.0001,0.0005,0.002'), ',', fixed=TRUE)[[1L]]),
  alpha = 0.05,
  folds = 5L,
  seed = 20261019L,
  atom_probability = 0.35,
  L = 0.40,
  U = 0.90,
  s = 1,
  sp2_at_75 = 31795402 / 38109275,
  truth_tolerance = 1e-8
)

v14_assert <- function(cond, msg) if (!isTRUE(cond)) stop(msg,call.=FALSE)
v14_clip <- function(x) pmin(1-1e-7,pmax(1e-7,x))
v14_expit <- function(x) stats::plogis(x)
v14_log <- function(...) cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '| V14 |', ..., '\n')

# The X distribution has a point mass at 30, plus a continuous component.
# Age may be interpreted as a constructed numeric covariate, not UCI ages.
v14_marginals <- function(x, direction='lower_zero') {
  t <- (x-30)/45
  v14_assert(all(is.finite(t) & t>=-1e-10 & t<=1+1e-10), 'x outside [30,75]')
  p <- 0.30 - 0.10*t
  se1 <- 0.80 + 0.05*t
  sp1 <- (32 - 4*t)/35
  se2 <- 0.80 - 0.05*t
  sp2 <- 32/35 + (v14_options$sp2_at_75-32/35)*t
  if (identical(direction,'upper_zero')) {
    tmp <- se1; se1<-se2; se2<-tmp
    tmp <- sp1; sp1<-sp2; sp2<-tmp
  }
  q1 <- p*se1+(1-p)*(1-sp1)
  q2 <- p*se2+(1-p)*(1-sp2)
  stopifnot(all(p>0 & p<1 & se1>0 & se1<1 & se2>0 & se2<1),
            all(sp1>0 & sp1<1 & sp2>0 & sp2<1))
  data.frame(x=x,t=t,p=p,se1=se1,sp1=sp1,se2=se2,sp2=sp2,q1=q1,q2=q2)
}

# Population integrals are computed using the known mixture distribution:
#   P_X f = 0.35 f(30) + 0.65 int_0^1 f(30+45 t) dt.
# The lower/upper truth uses the EXACT nonsmooth envelope definition.
v14_expectation <- function(fun) {
  a <- v14_options$atom_probability
  a*fun(0) + (1-a)*stats::integrate(fun, lower=0,upper=1,
     rel.tol=1e-12, abs.tol=1e-12, subdivisions=2000L)$value
}

v14_truth <- function(direction='lower_zero') {
  calc <- function(t) v14_marginals(30+45*t,direction)
  C1 <- v14_expectation(function(t){r<-calc(t);r$p*(1-r$p)*(r$se1+r$sp1-1)})
  C2 <- v14_expectation(function(t){r<-calc(t);r$p*(1-r$p)*(r$se2+r$sp2-1)})
  V1 <- v14_expectation(function(t){r<-calc(t);(1-r$p)*r$q1})
  V2 <- v14_expectation(function(t){r<-calc(t);(1-r$p)*r$q2})
  b0 <- C1*V2-C2*V1
  hh <- function(t){r<-calc(t);C1*(r$p-r$q2)-C2*(r$p-r$q1)}
  L <- v14_options$L; U <- v14_options$U
  lower <- b0+v14_expectation(function(t) pmin(L*hh(t),U*hh(t)))
  upper <- b0+v14_expectation(function(t) pmax(L*hh(t),U*hh(t)))
  minR <- function(j) {
    v14_expectation(function(t) {
       r<-calc(t);q<-if(j==1L)r$q1 else r$q2
       L0<-r$p*(1-q);U0<-(1-r$p)*q
       pmin(L*L0+(1-L)*U0, U*L0+(1-U)*U0)
    })
  }
  data.frame(direction=direction,C1=C1,C2=C2,V1=V1,V2=V2,
             b0=b0,truth_lower_B=lower,truth_upper_B=upper,
             lower_reference_test1=minR(1L),lower_reference_test2=minR(2L),
             mass_at_exact_h_tie=v14_options$atom_probability)
}

# Exact boundary certificate: derived symbolically by integrating polynomials.
# Let z=Sp2(75). Then B(c=L) = -13*(38109275*z-31795402)/20160000000.
# At z*=31795402/38109275, h(t)>=0 on t in [0,1], and P_X(h=0)=0.35.
# Hence lowerB=B(L)=0 exactly and upperB=556233219113/548773560000000.
# Interchanging test names negates/switches these two endpoints.
v14_truth_exact <- data.frame(
  direction=c('lower_zero','upper_zero'),
  truth_lower_B=c(0,-556233219113/548773560000000),
  truth_upper_B=c(556233219113/548773560000000,0),
  mass_at_exact_h_tie=c(0.35,0.35))

v14_generate <- function(n,rho,direction) {
  a <- v14_options$atom_probability
  atom <- stats::runif(n)<a
  x <- ifelse(atom,30,30+45*stats::runif(n))
  pa <- v14_marginals(x,direction)
  D <- stats::rbinom(n,1L,pa$p)
  r1 <- ifelse(D==1,pa$se1,1-pa$sp1)
  r2 <- ifelse(D==1,pa$se2,1-pa$sp2)
  u1 <- stats::runif(n)
  u2 <- ifelse(stats::runif(n)<rho,u1,stats::runif(n))
  data.frame(X=x,D=D,T1=as.integer(u1<r1),T2=as.integer(u2<r2),
             p_truth=pa$p,q1_truth=pa$q1,q2_truth=pa$q2)
}

# Replicates the brglm2 mean-bias-reduced out-of-fold logistic spline learner,
# without invoking the expensive V13.4 one-click research stages.
v14_nuisance <- function(dat, mode, seed) {
  if (identical(mode,'oracle')) {
    return(list(p=dat$p_truth,q1=dat$q1_truth,q2=dat$q2_truth))
  }
  v14_assert(identical(mode,'crossfit'), "Estimator must be 'oracle' or 'crossfit'")
  if(!requireNamespace('brglm2',quietly=TRUE))stop('Install brglm2 for crossfit mode')
  if(!requireNamespace('splines',quietly=TRUE))stop('splines unavailable')
  n <- nrow(dat); K<-v14_options$folds
  set.seed(seed)
  fold <- sample(rep(seq_len(K),length.out=n),size=n)
  pr <- matrix(NA_real_,n,3L)
  for(k in seq_len(K)){
    train<-dat[fold!=k,,drop=FALSE]
    val<-dat[fold==k,,drop=FALSE]
    for(j in seq_along(c('D','T1','T2'))){
      f<-stats::as.formula(paste0(c('D','T1','T2')[j],
                ' ~ splines::ns(X, df=3)'))
      fit<-stats::glm(f,data=train,family=stats::binomial(),
                method=brglm2::brglmFit,type='AS_mean',
                control=stats::glm.control(maxit=200))
      if(!isTRUE(fit$converged))stop('Nonconverged out-of-fold logistic model')
      pr[fold==k,j] <- v14_clip(stats::predict(fit,newdata=val,type='response'))
    }
  }
  if(any(!is.finite(pr)))stop('Invalid out-of-fold probability')
  list(p=pr[,1],q1=pr[,2],q2=pr[,3])
}

v14_soft <- function(h,tau,upper=FALSE){
  L<-v14_options$L; U<-v14_options$U
  if (!upper) {
    a <- -L*h/tau; b <- -U*h/tau
    mm <- pmax(a,b)
    return(-tau*(mm+log(exp(a-mm)+exp(b-mm))))
  }
  a <- L*h/tau; b <- U*h/tau
  mm <- pmax(a,b)
  tau*(mm+log(exp(a-mm)+exp(b-mm)))
}
v14_soft_grad <- function(h,tau,upper=FALSE){
  L<-v14_options$L; U<-v14_options$U
  if(!upper){
    pL<-stats::plogis((U-L)*h/tau)
    return(L*pL+U*(1-pL))
  }
  pU<-stats::plogis((U-L)*h/tau)
  L*(1-pU)+U*pU
}

# This is the ORIGINAL V13.4 one-step smoothed envelope/EIF algebra.
# The only distinction between modes is true vs cross-fitted p,q inputs.
v14_envelope_estimate <- function(dat,nu,tau,upper=FALSE) {
  n<-nrow(dat);D<-dat$D;T1<-dat$T1;T2<-dat$T2
  p<-nu$p;q1<-nu$q1;q2<-nu$q2
  eD<-D-p;e1<-T1-q1;e2<-T2-q2
  uC1<-eD*e1;uC2<-eD*e2
  C1<-mean(uC1);C2<-mean(uC2)
  v1<-(1-p)*q1-q1*eD+(1-p)*e1
  v2<-(1-p)*q2-q2*eD+(1-p)*e2
  V1<-mean(v1);V2<-mean(v2)
  phiC1<-uC1-C1;phiC2<-uC2-C2
  phiV1<-v1-V1;phiV2<-v2-V2
  h<-C1*(p-q2)-C2*(p-q1)
  b0<-C1*V2-C2*V1
  xi<-(C1-C2)*eD+C2*e1-C1*e2
  g<-v14_soft(h,tau,upper)
  w<-v14_soft_grad(h,tau,upper)
  alpha1<-mean(w*(p-q2));alpha2<-mean(w*(p-q1))
  theta<-b0+mean(g+w*xi)
  phi<-(V2+alpha1)*phiC1-(V1+alpha2)*phiC2+
       C1*phiV2-C2*phiV1+(g-mean(g))+w*xi
  se<-stats::sd(phi)/sqrt(n)
  crit<-stats::qnorm(1-v14_options$alpha/(2L*length(v14_options$taus)))
  bound<-if(!upper)theta-crit*se else theta+crit*se
  c(theta=theta,se=se,confidence_bound=bound,
    certified=as.integer(is.finite(bound) && if(!upper)bound>0 else bound<0))
}

v14_run <- function(){
  cfg<-v14_options
  stopifnot(cfg$estimator%in%c('oracle','crossfit'),
    length(cfg$reps)==1L,is.finite(cfg$reps),cfg$reps>=10,
    all(cfg$ns>=100),all(cfg$rhos>=0 & cfg$rhos<=1),
    all(cfg$taus>0),length(cfg$taus)>=1L)
  dir.create(cfg$out_dir,recursive=TRUE,showWarnings=FALSE)
  v14_assert(dir.exists(cfg$out_dir),'Output directory cannot be created')
  truths<-do.call(rbind,lapply(v14_truth_exact$direction,v14_truth))
  truths$exact_lower_B <- v14_truth_exact$truth_lower_B
  truths$exact_upper_B <- v14_truth_exact$truth_upper_B
  stopifnot(all(abs(truths$truth_lower_B-truths$exact_lower_B)<cfg$truth_tolerance),
            all(abs(truths$truth_upper_B-truths$exact_upper_B)<cfg$truth_tolerance),
            all(truths$lower_reference_test1>0 & truths$lower_reference_test2>0))
  utils::write.csv(truths,file.path(cfg$out_dir,'V14_EXACT_BOUNDARY_TRUTHS.csv'),row.names=FALSE)
  settings<-data.frame(key=names(cfg),value=vapply(cfg,function(z) paste(z,collapse=','),character(1)))
  utils::write.csv(settings,file.path(cfg$out_dir,paste0('V14_SETTINGS_',cfg$estimator,'.csv')),row.names=FALSE)
  allrows<-list(); itr<-0L
  for (direction in v14_truth_exact$direction) for (n in cfg$ns) for (rho in cfg$rhos) {
    v14_log('Start ',direction,' n=',n,' rho=',rho,' reps=',cfg$reps,' mode=',cfg$estimator)
    for (rep in seq_len(cfg$reps)) {
      # Deterministic separate streams per cell and replication.
      seed<-cfg$seed+as.integer(direction=='upper_zero')*100000000L+
          as.integer(n)*7001L+as.integer(round(rho*100))*100003L+rep
      set.seed(seed)
      obj<-tryCatch({
        dat<-v14_generate(n,rho,direction)
        nu<-v14_nuisance(dat,cfg$estimator,seed+334L)
        est<-do.call(rbind,lapply(cfg$taus,function(tau){
          lo<-v14_envelope_estimate(dat,nu,tau,upper=FALSE)
          hi<-v14_envelope_estimate(dat,nu,tau,upper=TRUE)
          data.frame(tau=tau,lower_theta=lo['theta'],lower_se=lo['se'],
            lower_LCB=lo['confidence_bound'],T1_certificate=lo['certified'],
            upper_theta=hi['theta'],upper_se=hi['se'],
            upper_UCB=hi['confidence_bound'],T2_certificate=hi['certified'])
        }))
        data.frame(direction=direction,n=n,rho=rho,rep=rep,
          estimation_mode=cfg$estimator,success=TRUE,failure_reason='',est)
      },error=function(e){
        data.frame(direction=direction,n=n,rho=rho,rep=rep,
          estimation_mode=cfg$estimator,success=FALSE,failure_reason=conditionMessage(e),
          tau=cfg$taus,lower_theta=NA_real_,lower_se=NA_real_,
          lower_LCB=NA_real_,T1_certificate=NA_integer_,
          upper_theta=NA_real_,upper_se=NA_real_,
          upper_UCB=NA_real_,T2_certificate=NA_integer_)
      })
      itr<-itr+1L;allrows[[itr]]<-obj
      if(rep%%100L==0L)v14_log('Progress ',direction,' n=',n,' rho=',rho,
                                      ' rep ',rep,'/',cfg$reps)
    }
    # Crash resilience: write every complete (direction,n,rho) cell.
    idx<-which(vapply(allrows,function(u)u$direction[1]==direction &&
        u$n[1]==n && abs(u$rho[1]-rho)<1e-12,logical(1)))
    z<-do.call(rbind,allrows[idx]); rownames(z)<-NULL
    filename<-sprintf('V14_RAW_%s_n%d_rho%03d_%s.csv',direction,n,
                     as.integer(round(rho*100)),cfg$estimator)
    utils::write.csv(z,file.path(cfg$out_dir,filename),row.names=FALSE)
  }
  raw<-do.call(rbind,allrows);rownames(raw)<-NULL
  summary_rows<-list();j<-0L
  for(direction in v14_truth_exact$direction)for(n in cfg$ns)for(rho in cfg$rhos)for(tau in cfg$taus){
    d<-raw[raw$direction==direction & raw$n==n & abs(raw$rho-rho)<1e-12 & abs(raw$tau-tau)<1e-12,,drop=FALSE]
    ok<-d$success
    wrong<-if(direction=='lower_zero')as.integer(d$T1_certificate==1) else as.integer(d$T2_certificate==1)
    anycert<-as.integer(d$T1_certificate==1 | d$T2_certificate==1)
    rate<-sum(wrong[ok],na.rm=TRUE)/nrow(d)
    j<-j+1L
    summary_rows[[j]]<-data.frame(direction=direction,n=n,rho=rho,tau=tau,
      estimator=cfg$estimator,replications=nrow(d),success=sum(ok),failures=sum(!ok),
      false_direction_certificates=sum(wrong[ok],na.rm=TRUE),
      false_direction_rate_all_reps=rate,
      any_certificate_rate_all_reps=sum(anycert[ok],na.rm=TRUE)/nrow(d),
      false_direction_mcse=sqrt(rate*(1-rate)/nrow(d)),
      truth_lower_B=v14_truth_exact$truth_lower_B[v14_truth_exact$direction==direction],
      truth_upper_B=v14_truth_exact$truth_upper_B[v14_truth_exact$direction==direction])
  }
  report<-do.call(rbind,summary_rows)
  utils::write.csv(report,file.path(cfg$out_dir,paste0('V14_BOUNDARY_SIZE_SUMMARY_',cfg$estimator,'.csv')),row.names=FALSE)
  writeLines(c('V14 BOUNDARY VALIDATION - PILOT DIAGNOSTIC',
     'Two analytically calibrated population boundary nulls with 35% mass at h(X)=0.',
     paste('Estimator mode:',cfg$estimator),
     paste('Monte Carlo replications per cell:',cfg$reps),
     'Each direction is a distinct boundary null; these are NOT 120 original V13 cells.',
     'Oracle mode isolates smoothing/test geometry; it does not certify fitted nuisance inference.',
     'Crossfit mode implements same 5-fold brglm2 logistic-spline recipe, but its consistency must be justified.',
     'Bonferroni correction uses alpha/(2*number of taus) for BOTH directions.',
     'Failure-adjusted fractions use all simulated replications as denominator.',
     'Fixed tau can be substantially conservative near the boundary.',
     'This is not a proof of uniform validity or a journal acceptance certificate.'),
     file.path(cfg$out_dir,paste0('V14_INTERPRETATION_AND_LIMITS_',cfg$estimator,'.txt')))
  v14_log('COMPLETED; files in ',cfg$out_dir)
  invisible(report)
}

v14_run()
