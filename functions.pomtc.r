###########################################################################
#       Model-based clustering using a mixture of transitional models
#       POM+transitional+interactions+time
#       Functions that implement MH sampler
#       Roy Costilla
#       Apr18
#       likelihood in C++ (theta.pomct and Zrc_tt)
##############################################################################


# Loading R packages
mypackages=c("MASS","coda","gtools","truncnorm","xtable")
mypackages.loaded=sapply(mypackages, require, character.only=T)
print(mypackages.loaded)
if (sum(mypackages.loaded)<length(mypackages)) {
  install.packages(mypackages[!mypackages.loaded], dependencies = T)
  sapply(mypackages[!mypackages.loaded], require, character.only=T)
}


# Main wrapper function to be passed to lapply
pomtc=function(init, mydata, fpar, qpar){
  cat("######################################### \n")
  cat("Init", init, "\n")
  param.start = Inits[[init]]
  state <- init.state(param.start , mydata, fpar, qpar)
  stopifnot(!any(is.na(state)))
  cat("mu=",round(state$params$mu,1), " alpha=",round(state$params$alpha,1)," sigma2.beta=",round(state$params$sigma2.beta,1)," pi=",round(state$params$pi,1), "  log.like=",round(state$log.like,1),"  log.prior=",round(state$log.prior,1),"\n", sep=' ')  
  print(round((state$params$gamma),1))  
  print(round(state$params$beta,1))
  chain <- run.chain(nburn=nburn, nthin=nthin, nstore=nstore,
                     state=state, data=mydata, fpar=fpar, qpar=qpar)
  cat("\n \n")
  chain
}


# Calculates theta given parameters
##### If you need to compile again (Rtools) uncomment next line
# system("R CMD SHLIB theta.pomtc.cpp --clean --preclean")
switch(Sys.info()[['sysname']],
       Windows= {dyn.load("theta.pomtc.dll")},
       Linux  = {dyn.load("theta.pomtc.so") })
if (is.loaded("theta_pomtc") ) cat("function theta_pomtc is loaded", '\n')	
theta.pomtc <- function(mu, alpha, beta, gamma,fpar) {
  q=fpar$q; p=fpar$p; R=fpar$R
  theta.temp <- array(1, dim=c(p,q,q,R))
  # Calling to C++ 
  out <- .C("theta_pomtc", mu=as.double(mu), alpha=as.double(alpha), beta=as.double(beta), gamma=as.double(gamma), R=as.integer(R), p=as.integer(p),  q=as.integer(q), thetacdf=as.double(theta.temp), theta=as.double(theta.temp))
  theta.temp <- array(out$theta, dim=c(R,q,q,p))
  stopifnot(theta.temp>=0) # theta is positive
  stopifnot(  (apply(theta.temp,c(1,2,4),sum)-1) <=1e-8) #theta adds to 1 over k (or close!)
  theta.temp
}

# Z for transitional model with column effects
##### If you need to compile again (Rtools) uncomment next line
#system("R CMD SHLIB Zrc_tt.cpp --clean --preclean")
switch(Sys.info()[['sysname']],
       Windows= {dyn.load("Zrc_tt.dll")},
       Linux  = {dyn.load("Zrc_tt.so") })
if (is.loaded("Zrc_tt") ) cat("function Zrc_tt is loaded", '\n')	

Z.rc.tt <- function(theta.temp, pi.temp ,fpar, data) {
  n=fpar$n
  p=fpar$p
  q=fpar$q
  R=fpar$R
  y=data$y
  # Incomplete information Log-likehood expressed as a function of Z 
  Z.temp=matrix(99,n,R) # setting up Z
  # C++ function that returns numerator of Z in logs
  out <- .C("Zrc_tt", as.double(theta.temp), as.integer(y), as.double(pi.temp), as.integer(n), as.integer(p), as.integer(q), as.integer(R), Z=as.double(Z.temp))
  Z.temp = t(matrix(out$Z, R, n)) # transposing C output vector (R rows x n cols)  
  stopifnot(!is.na(Z.temp))
  Z.temp
}

# Incomplete information log-like (different from EM one, where is the complete!!)
log.like.func <- function(state, data , fpar) {
  n=fpar$n
  p=fpar$p
  q=fpar$q
  R=fpar$R
  y=data$y
  theta = theta.pomtc(state$params$mu,state$params$alpha,state$params$beta,state$params$gamma,fpar)
  pi=state$params$pi
  # Incomplete information Log-likehood expressed as a function of Z 
  Z=matrix(99,n,R) # setting up Z
  # C++ function that returns numerator of Z in logs
  out <- .C("Zrc_tt", as.double(theta), as.integer(y), as.double(pi), as.integer(n), as.integer(p), as.integer(q), as.integer(R), Z=as.double(Z))
  Z = t(matrix(out$Z, R, n)) # transposing C output vector (R rows x n cols)  
  stopifnot(!is.na(Z))
  log.like <- sum(log(apply(exp(Z),1,sum)))
  log.like
}

# Initializing MCMC state
init.state <- function(param, data, fpar, qpar) {
  stopifnot(fpar$R>=1)
  stopifnot(sum(param$pi)-1<= 1e-10)
  state <- list(params=param,
                log.like=NA, log.prior=NA,
                accepted=NA)
  state$log.prior <- log.prior.func(state, fpar)
  state$log.like <- log.like.func(state, data, fpar)  
  state$accepted <- rep(TRUE,length(fpar$vars))
  names(state$accepted) <-fpar$vars
  return(state)
}

# Priors
log.prior.func <- function(state, fpar) {
  # Uniform doubly-truncated priors for cut points  mu ~ U[-20,20]I[mu_k-1,mu_k+1]
  # k=1
  log.prior.mu= dunif(state$params$mu[1], min=0, max=state$params$mu[2], log = T)
  # 2 <= k <= q-2
  if (q>3) for (k in 2:(q-2)) log.prior.mu = log.prior.mu +
      dunif(state$params$mu[k], min=state$params$mu[k-1], max=state$params$mu[k], log = T)
  # k= q-1
  log.prior.mu = dunif(state$params$mu[q-1], min=state$params$mu[q-2], max=20, log = T) +
    log.prior.mu
  log.prior.sigma2.mu = dgamma(1/state$params$sigma2.mu, shape=fpar$a, rate=fpar$b, log=TRUE)
  
  # degenerate normal prior for beta_rk' 
  scalingf=(q-1)/q
  if (fpar$R==1) log.prior.beta = sum(dnorm(state$params$beta[-q],
                                            log=T,sd=scalingf*sqrt(state$params$sigma2.beta))) else 
                                              log.prior.beta = sum(apply(state$params$beta[,-q],2,dnorm,mean=0,log=T,sd=scalingf*sqrt(state$params$sigma2.beta)))
  log.prior.sigma2.beta = sum(dgamma(1/state$params$sigma2.beta, shape=fpar$a, rate=fpar$b, log=TRUE))
  
  # gamma
  # degenerate normal prior
  scalingf=(p-2)/(p-1)
  log.prior.gamma = sum(dnorm(state$params$gamma[c(-1,-2)], mean=0, 
                              sd=scalingf*sqrt(state$params$sigma2.gamma), log=TRUE))
  
  log.prior.sigma2.gamma = dgamma(1/state$params$sigma2.gamma, shape=fpar$a, rate=fpar$b, log=TRUE)
  
  #alpha and pi
  if (fpar$R==1) {
    log.prior.alpha=0
    log.prior.sigma2.alpha=0
    log.prior.pi=0
  } else {
    log.prior.alpha = sum(dnorm(state$params$alpha, mean=0, sd=sqrt(state$params$sigma2.alpha), log=TRUE))
    log.prior.sigma2.alpha = dgamma(1/state$params$sigma2.alpha, shape=fpar$a, rate=1/fpar$b, log=TRUE)
    log.prior.pi = log(ddirichlet(state$params$pi, fpar$phi))
  }
  
  # Adding up all prior contributions
  log.prior.func =  log.prior.mu+log.prior.sigma2.mu+
    log.prior.alpha+log.prior.sigma2.alpha+
    log.prior.beta+log.prior.sigma2.beta+
    log.prior.gamma+log.prior.sigma2.gamma+
    log.prior.pi
  log.prior.func
}

update.chain <- function(state, data, fpar, qpar, par) {
  n=fpar$n
  p=fpar$p
  q=fpar$q
  R=fpar$R
  old.state <- state
  ############################      Sigma2's
  # Log-Normal RW proposal. Includes Jacobian for r
  if ( grepl("sigma2",par) ) {   #if "sigma2" is part of string "par"
    if ( grepl("beta",par)  ) {
      state$params$sigma2.beta <-  rlnorm(fpar$R, mean=log(state$params$sigma2.beta), sdlog=qpar$proposal.var.beta)
      log.proposal.ratio <- sum(log(old.state$params$sigma2.beta) - log(state$params$sigma2.beta))
    } else{
      if ( grepl("gamma",par)  ) {
        state$params$sigma2.gamma <-  rlnorm(1, mean=log(state$params$sigma2.gamma), sdlog=qpar$proposal.var.gamma)
        log.proposal.ratio <- log(old.state$params$sigma2.gamma) - log(state$params$sigma2.gamma)
      } else {
        if( grepl("alpha",par) ) {
          if (R>1) {  
            state$params$sigma2.alpha <-  rlnorm(1, mean=log(state$params$sigma2.alpha), sdlog=qpar$proposal.var.alpha)
            log.proposal.ratio <- log(old.state$params$sigma2.alpha) - log(state$params$sigma2.alpha)
          } else {
            log.proposal.ratio <- 0
          }
        }else {  #mu
          state$params$sigma2.mu <-  rlnorm(1, mean=log(state$params$sigma2.mu), sdlog=qpar$proposal.var.mu)
          log.proposal.ratio <- log(old.state$params$sigma2.mu) - log(state$params$sigma2.mu)
        }
      }
    }
    state$log.like <- old.state$log.like  # likelihood only depends on beta and mu
    state$log.prior <- log.prior.func(state, fpar) #flat=TRUE
    # compute log acceptance ratio
    log.like.ratio <- state$log.like - old.state$log.like
    log.prior.ratio <- state$log.prior - old.state$log.prior
    log.r <- log.like.ratio + log.prior.ratio + log.proposal.ratio
    # accept or reject the state
    u = runif(1)      
    if(u <exp(log.r)) {
      # accept
      state$accepted[par] <- TRUE
    } else {
      # reject
      state <-old.state
      state$accepted[par] <- FALSE
    }
  } else {
    ############################      beta 
    if (par=='beta') {
      # beta_rk', single component update
      r=sample(1:R, 1)
      kdash= sample(1:(q-1), 1)
      state$params$beta[r,kdash]=rnorm(1, mean=state$params$beta[r,kdash], sd=qpar$proposal.beta)
      # Last beta is the sum of the others 
      if (fpar$R==1) state$params$beta[q] = sum(state$params$beta[-q]) else 
        state$params$beta[,q]=-apply(state$params$beta[,-q],1,sum)
      
      # likelihood and prior for new state
      state$log.like <- log.like.func(state, data, fpar) #ccode=TRUE  check this!!
      state$log.prior <- log.prior.func(state, fpar)
      # compute log acceptance ratio
      log.like.ratio <- state$log.like - old.state$log.like
      log.prior.ratio <- state$log.prior - old.state$log.prior
      #log.proposal.ratio <- log(state$params$tau2) - log(old.state$params$tau2)
      log.proposal.ratio <- 0
      log.r <- log.like.ratio + log.prior.ratio + log.proposal.ratio
      # accept or reject the state
      u = runif(1)      
      if(u <exp(log.r)) {
        # accept
        state$accepted[par] <- TRUE
      } else {
        # reject
        state <-old.state
        state$accepted[par] <- FALSE
      }
    } else {
      ############################      gamma 
      if (par=='gamma') {
        
        # p/2 components update
        j= sample(1:(p-2), floor(p/2))+2
        state$params$gamma[j]=rnorm(length(j), mean=state$params$gamma[j], sd=qpar$proposal.gamma)
        state$params$gamma[2]= -sum(state$params$gamma[3:p])
        
        # likelihood and prior for new state
        state$log.like <- log.like.func(state, data, fpar) 
        state$log.prior <- log.prior.func(state, fpar)
        
        # compute log acceptance ratio
        log.like.ratio <- state$log.like - old.state$log.like
        log.prior.ratio <- state$log.prior - old.state$log.prior
        log.proposal.ratio <- 0
        log.r <- log.like.ratio + log.prior.ratio + log.proposal.ratio
        
        # accept or reject the state
        u = runif(1)      
        if(u <exp(log.r)) {
          # accept
          state$accepted[par] <- TRUE
        } else {
          # reject
          state <-old.state
          state$accepted[par] <- FALSE
        }
      }  else {
        ############################      Alpha
        if (par=='alpha') {
          # update alpha only for R>1, otherwise do nothing (alpha=0, no clustering)
          if (R>1) {
            # single component update
            g= sample(1:(R), 1)
            alpha.star = rnorm(1, mean=state$params$alpha[g], sd=qpar$proposal.alpha)
            state$params$alpha[g]=alpha.star
            # likelihood and prior for new state
            state$log.like <- log.like.func(state, data, fpar) #ccode=TRUE  check this!!
            # compute log acceptance ratio
            log.like.ratio <- state$log.like - old.state$log.like
            log.prior.ratio <- state$log.prior - old.state$log.prior
            log.proposal.ratio <- 0
            log.r <- log.like.ratio + log.prior.ratio + log.proposal.ratio
            # accept or reject the state
            u = runif(1)  
            if(u <exp(log.r)) {
              # accept
              state$accepted[par] <- TRUE
            } else {
              # reject
              state <-old.state
              state$accepted[par] <- FALSE
            }
          }
        } else {
          ############################      Mu  (block)
          # mu is sample in block from truncated uniform
          if (par=='mu') {
            tau=qpar$proposal.mu
            mu=c(-Inf,state$params$mu,Inf) 
            k=ifelse(q==3, 2,sample((2:(q-1)),1))
            index=k+1 # index in R starts on 1
            mu.star=runif(1,max(mu[index]-tau,mu[index-1]), min(mu[index]+tau,mu[index+1]))
            state$params$mu[k] <- mu.star
            
            # New log-like and priors
            state$log.prior <- log.prior.func(state, fpar)
            state$log.like <- log.like.func(state, data, fpar)  
            
            # compute log acceptance ratio
            log.like.ratio <- state$log.like - old.state$log.like
            log.prior.ratio <- state$log.prior - old.state$log.prior
            log.proposal.ratio <- -log(min(mu[index]+tau,mu[index+1]) - max(mu[index]-tau,mu[index-1]))
            +log(min(mu.star+tau,mu[index+1]) - max(mu.star-tau,mu[index-1]))
            log.r <- log.like.ratio + log.prior.ratio + log.proposal.ratio
            # accept or reject the state
            u = runif(1)  
            if(u <exp(log.r)) {
              # accept
              state$accepted[par] <- TRUE
            } else {
              # reject
              state <-old.state
              state$accepted[par] <- FALSE
            }
          } else {
            ######################        pi  (2 at a time for any R)
            # update alpha only for R>1, otherwise do nothing (alpha=0)
            if (par=='pi') {
              if (R>1) {
                # pick 2 at a time
                pi=old.state$params$pi
                pi.star=old.state$params$pi
                r=sort(sample(1:fpar$R,2,replace=FALSE))
                
                # Proportion of total pi to be change pi[r[1]]+pi[r[2]]
                w=pi[r[1]]/(pi[r[1]]+pi[r[2]])
                
                # Big enough to get negative values for logit (pi.star is less pi[r]
                logit.w.star = rnorm(1, log(w)-log(1-w) ,qpar$proposal.sd.pi)
                w.star = 1/(1+exp(-logit.w.star))
                pi.star[r[1]] = w.star*(pi[r[1]]+pi[r[2]])
                pi.star[r[2]] = (1-w.star)*(pi[r[1]]+pi[r[2]])
                state$params$pi[r[1]] = pi.star[r[1]]
                state$params$pi[r[2]] = pi.star[r[2]]
                
                # New log-like and priors
                state$log.like <- log.like.func(state, data, fpar) 
                state$log.prior <- log.prior.func(state, fpar)
                
                # compute log acceptance ratio
                log.like.ratio <- state$log.like - old.state$log.like
                log.prior.ratio <- state$log.prior - old.state$log.prior
                log.proposal.ratio <- ( log(w.star)+log(1-w.star) )
                -( log(w)+log(1-w) )
                log.r <- log.like.ratio + log.prior.ratio + log.proposal.ratio
                
                # accept or reject the state
                u = runif(1)  
                if(u <exp(log.r)) {
                  # accept
                  state$accepted[par] <- TRUE
                  #cat("accepted, pi=", round(state$params$pi,2), '\n')
                } else {
                  # reject
                  state <-old.state
                  state$accepted[par] <- FALSE
                }
              }
            }
          }
          
        }
      }
    }
  }
  state
}

# Wrapper function to run MCMC chain
run.chain <- function(nburn, nthin, nstore,
                      state, data, fpar, qpar) {
  n=fpar$n
  p=fpar$p
  q=fpar$q
  R=fpar$R
  cat('Burn-in ', nburn, '\n')
  for(i in 1:(nburn/9)) {
    state <- update.chain(state, data, fpar, qpar, par="sigma2.mu")
    state <- update.chain(state, data, fpar, qpar, par="mu")
    state <- update.chain(state, data, fpar, qpar, par="sigma2.alpha")
    state <- update.chain(state, data, fpar, qpar, par="alpha")
    state <- update.chain(state, data, fpar, qpar, par="sigma2.beta")
    state <- update.chain(state, data, fpar, qpar, par="beta")
    state <- update.chain(state, data, fpar, qpar, par="gamma")
    state <- update.chain(state, data, fpar, qpar, par="sigma2.gamma")
    state <- update.chain(state, data, fpar, qpar, par="pi")
    if( i/1000==floor(i/1000) )  cat('+')
  }
  cat('+ \n')
  cat("State after burn-in=", '\n')
  cat("mu=",round(state$params$mu,1), " alpha=",round(state$params$alpha,1)," sigma2.beta=",round(state$params$sigma2.beta,1)," pi=",round(state$params$pi,2), "  log.like=",round(state$log.like,1),"  log.prior=",round(state$log.prior,1),"\n", sep=' ')
  cat('gamma=',round(state$params$gamma,1), '\n')
  cat('beta')
  print(round(state$params$beta,1))
  
  cat('\n')
  cat("Sampling=", nstore, 'thin=', nthin, 'total=', nstore*nthin , '\n')
  
  # chain to be stored
  npars = fpar$nfollow - length(fpar$vars)
  chain <- array(NA,dim=c(nstore,npars))
  # names for chain columns
  allnames <- c(paste("mu",1:(q-1),sep=''),"sigma2.mu")
  allnames <- c(allnames, paste("alpha",1:(R),sep=''), "sigma2.alpha")
  betas=paste("beta",1,'.',1:(q),sep='')
  if (R>1) for (r in 2:R) betas=c(betas,paste("beta",r,'.',1:(q),sep=''))
  allnames <- c(allnames, betas , paste("sigma2.beta",1:(R),sep=''))
  
  allnames <- c(allnames, paste("gamma",1:(p),sep=''),"sigma2.gamma")
  allnames <- c(allnames, paste("pi",1:(R),sep=''))
  allnames <- c(allnames, "llike","lpost")
  dimnames(chain)[[2]] <-allnames
  
  # Thinning (only if nthin>1, update states "nthin-1" times)
  #Update and save the "thin-th" iteration
  for(i in 1:(nstore/9)) {
    if (nthin>1) {
      for (t in 1:(nthin-1) ) {
        state <- update.chain(state, data, fpar, qpar, par="sigma2.mu")
        state <- update.chain(state, data, fpar, qpar, par="mu")
        state <- update.chain(state, data, fpar, qpar, par="sigma2.alpha")
        state <- update.chain(state, data, fpar, qpar, par="alpha")
        state <- update.chain(state, data, fpar, qpar, par="sigma2.beta")
        state <- update.chain(state, data, fpar, qpar, par="beta")
        state <- update.chain(state, data, fpar, qpar, par="sigma2.gamma")
        state <- update.chain(state, data, fpar, qpar, par="gamma")
        state <- update.chain(state, data, fpar, qpar, par="pi")
      }
    }
    
    # Update mu
    index=9*(i-1)+1
    for (thisvar in vars) {
      state <- update.chain(state, data, fpar, qpar, par=thisvar)
      chain[index,] <- c(state$params$mu, state$params$sigma2.mu,
                         state$params$alpha,state$params$sigma2.alpha,
                         as.vector(t(state$params$beta)),state$params$sigma2.beta,
                         state$params$gamma,state$params$sigma2.gamma,
                         state$params$pi,
                         state$log.like,state$log.like+state$log.prior)
      #,state$accept)
      index=index+1
    }
    # Printing symbols to monitor chain every 100 and 1000 iterations
    if(i/100==floor(i/100)) cat('.')
    if(i/1000==floor(i/1000)) {
      cat('\n', 'iteration=', i, '\n')
      cat("mu=",round(state$params$mu,1), " alpha=",round(state$params$alpha,1)," sigma2.beta=",round(state$params$sigma2.beta,1)," pi=",round(state$params$pi,2), "  log.like=",round(state$log.like,1),"  log.prior=",round(state$log.prior,1),"\n", sep=' ')
      cat('gamma=',round(state$params$gamma,1), '\n')
      cat('beta')
      print(round(state$params$beta,1))
    }
    
  }
  cat('\n')
  cat("mu=",round(state$params$mu,1), " alpha=",round(state$params$alpha,1)," sigma2.beta=",round(state$params$sigma2.beta,1)," pi=",round(state$params$pi,2), "  log.like=",round(state$log.like,1),"  log.prior=",round(state$log.prior,1),"\n", sep=' ')
  cat('gamma=',round(state$params$gamma,1), '\n')
  cat('beta')
  print(round(state$params$beta,1))
  return(chain)
}


# Entropy calculation
entropy <- function(this.Z) {
  this.Z[this.Z==0]=1e-6 # log(0) is -Inf 
  e= sum(apply(this.Z*log(this.Z),1,sum))
  e
}

# Cluster membership
group <- function(Z) {
  apply(Z,1,function(m) which(m==max(m)))
}
