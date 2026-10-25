# =============================================================================
# JBS V15 — INDEPENDENT NUISANCE-RATE DIAGNOSTIC FOR V14 BOUNDARY DGM
# Base research code, standalone. DOES NOT MODIFY V13 / V14 RESULTS.
# Goal: estimate population L2 errors of p(X),q1(X),q2(X) from models trained
# on approximately 4n/5 observations, matching a five-fold learner's fold size.
# The known V14 mixed covariate distribution permits near-exact deterministic
# numerical L2 integration (atom at age 30 plus equally spaced uniform grid).
# Empirical evidence cannot establish an asymptotic rate theorem on its own.
# =============================================================================
# In RStudio:
# install.packages('brglm2') # only if needed
# Sys.setenv(JBS_V15_REPS='50', JBS_V15_NS='250,500,858,1500,3000',
#            JBS_V15_DIR='C:/Users/user/Downloads/JBS_V15_NUISANCE_RESULTS')
# source(file.choose(), encoding='UTF-8') # select THIS .R file
# For a larger confirmatory run: use NEW output directory, JBS_V15_REPS='200'.
# =============================================================================

v15_cfg <- list(
  out_dir = Sys.getenv('JBS_V15_DIR', unset=file.path(getwd(),'JBS_V15_NUISANCE_RESULTS')),
  reps = as.integer(Sys.getenv('JBS_V15_REPS',unset='50')),
  ns = as.integer(strsplit(Sys.getenv('JBS_V15_NS',unset='250,500,858,1500,3000'),',',fixed=TRUE)[[1]]),
  models = c('linear','spline_df3','spline_df6'),
  directions = c('lower_zero','upper_zero'),
  folds = 5L,
  eval_grid = as.integer(Sys.getenv('JBS_V15_GRID',unset='2001')),
  atom_mass = 0.35,
  base_seed = 20261019L,
  sp2_end = 31795402 / 38109275,
  clip = 1e-7
)

v15_assert <- function(b,msg) if(!isTRUE(b)) stop(msg,call.=FALSE)
v15_log <- function(...) cat(format(Sys.time(),'%H:%M:%S'),' | V15 | ',..., '\n',sep='')
v15_expit <- function(z) stats::plogis(z)
v15_clip01 <- function(p) pmin(1-v15_cfg$clip,pmax(v15_cfg$clip,p))

# DGM copied from V14. The test dependence parameter rho is irrelevant for
# marginal p/q learner-rate diagnostics; we use rho=.25 to generate paired Y.
v15_truth <- function(x,direction) {
  tt <- (x-30)/45
  v15_assert(all(is.finite(tt) & tt>=-1e-9 & tt<=1+1e-9),'x out of DGM range')
  p <- 0.30-0.10*tt
  se1 <- 0.80+0.05*tt
  se2 <- 0.80-0.05*tt
  sp1 <- (32-4*tt)/35
  sp2 <- 32/35+(v15_cfg$sp2_end-32/35)*tt
  if(identical(direction,'upper_zero')){
    temp<-se1;se1<-se2;se2<-temp
    temp<-sp1;sp1<-sp2;sp2<-temp
  }
  q1 <- p*se1+(1-p)*(1-sp1)
  q2 <- p*se2+(1-p)*(1-sp2)
  data.frame(X=x,p=p,q1=q1,q2=q2,se1=se1,se2=se2,sp1=sp1,sp2=sp2)
}

v15_generate <- function(n,direction,rho=.25) {
  atom <- stats::runif(n) < v15_cfg$atom_mass
  x <- ifelse(atom,30,30+45*stats::runif(n))
  tt <- v15_truth(x,direction)
  D <- stats::rbinom(n,1L,tt$p)
  r1 <- ifelse(D==1,tt$se1,1-tt$sp1)
  r2 <- ifelse(D==1,tt$se2,1-tt$sp2)
  u1 <- stats::runif(n)
  u2 <- ifelse(stats::runif(n)<rho,u1,stats::runif(n))
  data.frame(X=x,D=D,T1=as.integer(u1<r1),T2=as.integer(u2<r2))
}

v15_fit_and_predict <- function(train,grid,model) {
  # brglm2 mean-bias reduction and natural splines exactly match the primary
  # nuisance fitting strategy; df6 and linear are sensitivity comparators.
  names_out <- c('D','T1','T2')
  preds <- matrix(NA_real_,nrow=nrow(grid),ncol=3L,
                  dimnames=list(NULL,c('p','q1','q2')))
  for(j in seq_along(names_out)){
    predictor <- switch(model,
      linear='X',
      spline_df3='splines::ns(X, df = 3)',
      spline_df6='splines::ns(X, df = 6)',
      stop('Unknown model'))
    f <- stats::as.formula(paste(names_out[j],'~',predictor))
    fit <- stats::glm(f,data=train,family=stats::binomial(),
                      method=brglm2::brglmFit,type='AS_mean',
                      control=stats::glm.control(maxit=200L))
    if(!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit))))
      stop(paste0('Nonconvergence in ',names_out[j],' / ',model))
    pp <- stats::predict(fit,newdata=grid,type='response')
    if(any(!is.finite(pp)))stop('Nonfinite nuisance predictions')
    preds[,j] <- v15_clip01(pp)
  }
  preds
}

v15_error_metrics <- function(pred,true) {
  # Grid row 1 = atom X=30, rows 2..(M+1) = Uniform[30,75] grid.
  a <- v15_cfg$atom_mass
  err2 <- function(j){
    sq <- (pred[,j]-true[[c('p','q1','q2')[j]]])^2
    a*sq[1]+(1-a)*mean(sq[-1])
  }
  ee <- sqrt(vapply(1:3,err2,numeric(1)))
  names(ee) <- c('ep','eq1','eq2')
  ee
}

v15_run <- function() {
  cfg <- v15_cfg
  v15_assert(requireNamespace('brglm2',quietly=TRUE),
             "Install 'brglm2' before running V15")
  v15_assert(requireNamespace('splines',quietly=TRUE),"splines unavailable")
  v15_assert(length(cfg$reps)==1L && is.finite(cfg$reps) && cfg$reps>=10,
             'JBS_V15_REPS must be an integer >=10')
  v15_assert(length(cfg$ns)>0 && all(is.finite(cfg$ns)) && all(cfg$ns>=100),
             'JBS_V15_NS must contain n >=100')
  v15_assert(cfg$eval_grid>=501,'Use integration grid >=501')
  dir.create(cfg$out_dir,recursive=TRUE,showWarnings=FALSE)
  v15_assert(dir.exists(cfg$out_dir),'Could not create output directory')
  out_raw <- file.path(cfg$out_dir,'V15_NUISANCE_RATES_RAW.csv')
  if(file.exists(out_raw))stop('Output RAW already exists. Choose NEW JBS_V15_DIR; no overwriting.')
  settings <- data.frame(setting=names(cfg),
                         value=vapply(cfg,function(z)paste(z,collapse=','),character(1)))
  utils::write.csv(settings,file.path(cfg$out_dir,'V15_SETTINGS.csv'),row.names=FALSE)
  x_eval <- c(30, 30+45*seq(0,1,length.out=cfg$eval_grid))
  eval_dat <- data.frame(X=x_eval)
  rows <- list();idx <- 0L
  for(direction in cfg$directions){
    truth <- v15_truth(x_eval,direction)
    for(n in cfg$ns){
      n_train <- floor(n*(cfg$folds-1)/cfg$folds)
      for(rep in seq_len(cfg$reps)){
        # A new sample per (direction,n,rep); all models evaluated on same data.
        seed <- cfg$base_seed + as.integer(direction=='upper_zero')*100000000L +
                n*2001L + rep
        set.seed(seed)
        train <- v15_generate(n_train,direction)
        for(model in cfg$models){
          obj <- tryCatch({
            pr <- v15_fit_and_predict(train,eval_dat,model)
            ee <- v15_error_metrics(pr,truth)
            ee
          },error=function(e)e)
          idx <- idx+1L
          if(inherits(obj,'error')){
            rows[[idx]] <- data.frame(direction=direction,n=n,train_n=n_train,
               rep=rep,model=model,success=FALSE,failure_reason=conditionMessage(obj),
               ep=NA_real_,eq1=NA_real_,eq2=NA_real_,
               prod1=NA_real_,prod2=NA_real_,sqrt_n_prod1=NA_real_,
               sqrt_n_prod2=NA_real_,sqrt_n_max_prod=NA_real_)
          }else{
            ep <- unname(obj['ep']); eq1<-unname(obj['eq1']);eq2<-unname(obj['eq2'])
            rows[[idx]] <- data.frame(direction=direction,n=n,train_n=n_train,
               rep=rep,model=model,success=TRUE,failure_reason='',
               ep=ep,eq1=eq1,eq2=eq2,
               prod1=ep*eq1,prod2=ep*eq2,sqrt_n_prod1=sqrt(n)*ep*eq1,
               sqrt_n_prod2=sqrt(n)*ep*eq2,
               sqrt_n_max_prod=sqrt(n)*ep*max(eq1,eq2))
          }
        }
        if(rep%%10L==0L){
          v15_log('direction=',direction,' n=',n,' rep=',rep,'/',cfg$reps)
          # Checkpoint to a separate file; the final RAW is not written until complete.
          if(rep%%20L==0L) utils::write.csv(do.call(rbind,rows),
              file.path(cfg$out_dir,'V15_CHECKPOINT_PARTIAL.csv'),row.names=FALSE)
        }
      }
    }
  }
  raw <- do.call(rbind,rows)
  utils::write.csv(raw,out_raw,row.names=FALSE)
  metrics <- c('ep','eq1','eq2','prod1','prod2','sqrt_n_prod1',
               'sqrt_n_prod2','sqrt_n_max_prod')
  groups <- split(raw,interaction(raw$direction,raw$n,raw$model,drop=TRUE))
  summary <- do.call(rbind,lapply(groups,function(g){
    good <- g[g$success,,drop=FALSE]
    ans <- data.frame(direction=g$direction[1],n=g$n[1],model=g$model[1],
          attempts=nrow(g),successful=nrow(good),failed=nrow(g)-nrow(good),
          success_rate=nrow(good)/nrow(g))
    for(v in metrics){
      z<-good[[v]]
      ans[[paste0('mean_',v)]] <- if(length(z)) mean(z) else NA_real_
      ans[[paste0('median_',v)]] <- if(length(z)) stats::median(z) else NA_real_
      ans[[paste0('sd_',v)]] <- if(length(z)>1) stats::sd(z) else NA_real_
    }
    ans
  }))
  rownames(summary)<-NULL
  summary <- summary[order(summary$direction,summary$model,summary$n),]
  utils::write.csv(summary,file.path(cfg$out_dir,'V15_NUISANCE_RATES_SUMMARY.csv'),row.names=FALSE)
  trend_groups <- split(summary,interaction(summary$direction,summary$model,drop=TRUE))
  trends <- do.call(rbind,lapply(trend_groups,function(g){
    out <- data.frame(direction=g$direction[1],model=g$model[1],
                      n_sizes=nrow(g),n_min=min(g$n),n_max=max(g$n))
    for(v in c('median_ep','median_eq1','median_eq2','median_sqrt_n_max_prod')){
      y<-g[[v]]; keep<-is.finite(y)&y>0
      out[[paste0('loglog_slope_',v)]] <- if(sum(keep)>=3)
        unname(stats::coef(stats::lm(log(y[keep])~log(g$n[keep])))[2]) else NA_real_
    }
    out
  }))
  rownames(trends)<-NULL
  utils::write.csv(trends,file.path(cfg$out_dir,'V15_RATE_TRENDS_EXPLORATORY.csv'),row.names=FALSE)
  writeLines(c(
    'V15 NUISANCE-RATE DIAGNOSTIC — RESEARCH EVIDENCE ONLY',
    'Expected target condition: sqrt(n)*||p_hat-p||_2*max_j||qj_hat-qj||_2 -> 0.',
    'True p,q from original V14 boundary DGM; X has 35% mass at 30.',
    'Models train on floor(4n/5) observations, as in one 5-fold training split.',
    'L2(P_X) is evaluated on an independent deterministic age mixture grid.',
    'These are separate holdout-training fits, NOT literal V14 five-fold predictions.',
    'The 3-df model is V13/V14 primary learner; linear/6-df are sensitivity checks.',
    'Fixed-dimensional misspecification can cause a nonvanishing approximation bias.',
    'Finite-sample trends do NOT prove asymptotic rates or uniform validity.',
    'Assess model failures and whether median sqrt(n)*max(product) decreases with n.',
    'Use output exclusively in a new V15 directory, preserving V13/V14 data.'
  ),file.path(cfg$out_dir,'V15_INTERPRETATION_README.txt'))
  part<-file.path(cfg$out_dir,'V15_CHECKPOINT_PARTIAL.csv')
  if(file.exists(part))file.remove(part)
  v15_log('COMPLETE: ',nrow(raw),' model-replication rows; ',nrow(summary),
          ' summary rows; output=',cfg$out_dir)
  invisible(list(raw=raw,summary=summary,trends=trends))
}

invisible(v15_run())
