#############################################################################
#       Model-based clustering using a mixture of transitional models
#       POM+transitional+interactions+time
#       Main program: 
#       (1) Simulate data 
#       (2) Run model
#       (3) Output postprocessing (relabel MCMC chains, etc )
#       Roy Costilla
#       Apr18
#############################################################################

#############################################################################
######################       (1) Simulate data        #######################

# Setting-up random number generator
rm(list = ls())
myseed=654321
set.seed(myseed)

# Loading R libraries and functions for model
source("functions.pomtc.r")

# Setting main parameters for simulation
thisdata='sim'
if (thisdata=='sim') {
  n=1000; p=15; q=5; G=R=3
  # Paramaters for Inverse Gamma prior
  mya=3;myb=1
}

# MCMC parameters
ninits=3
nburn=3*9e3
nstore=3*9e2
nthin=200

# Relabelling parameters
tol=1e-2
relabel=1

# Continue only if ordinal data (at least three categories)
stopifnot(q>=3)

# Setting true values for model paramters  
pi.true=rep(1/G,G)
if (G==3) alpha.true <- c(0,-1,1)
# Same probability for each ordinal response (1/q)
mu.true=qlogis((1:(q-1))/q, location=0, scale=1) 
cat('mu.true=', round(mu.true,2), 'alpha.true=', round(alpha.true,2), '\n \n')

# mu1=0 and no constraint in alpha
alpha.true=alpha.true-mu.true[1]
mu.true=mu.true-mu.true[1]

# since mu_k-alpha_r 
cat('mu.true=', round(mu.true,2), 'alpha.true=', round(alpha.true,2), '\n')

#beta
#sigma2.beta.true=0.05*(1:G)^2 # c(0.05, 0.2, 0.45)
sigma2.beta.true=c(0.2, 0.05, 0.45)
beta.true=matrix(-99,nrow=G,ncol=q)
for (r in 1:G) beta.true[r,] = sort(rnorm(q,mean=0,sd=sqrt(sigma2.beta.true[r])))
beta.true = beta.true-beta.true[,q]
# sum to zero constraint
beta.true[,q]=-apply(beta.true,1,sum)
cat('sigma2.beta.true=', round( sigma2.beta.true,2), ' beta \n')
print(round(beta.true,2))

# gamma
sigma2.gamma.true=2 
gamma.true = c(0,0, rnorm(p-2,mean=0,sd=sqrt(sigma2.gamma.true)))
gamma.true[2] = -sum(gamma.true)
cat('gamma.true=', round(gamma.true,2), '\n')

# Simulating y
r.true=sample(1:G, n, replace=T, prob=pi.true)  # distribution of cluster membership is pi.true
y=array(-99, dim=c(n,p+1))

# Beginning of time (assumption: all ordinal levels are equally likely)
y[,1]=matrix(sample(1:q,n, replace=T), ncol=1)
table(y[,1])
round(table(y[,1])/sum(table(y[,1])),2); mean(y[,1])

# Ocassion 1 onwards (2 to p+1)
fpar <- list(a=mya, b=myb, phi=3/2*rep(1,G), n=n, p=p ,q=q, R=G)
theta.pomtc.true=theta.pomtc(mu.true, alpha.true, beta.true, gamma.true, fpar)

# Transition probabilities averaged over time
cat('Cluster Transition probabilities averaged over time G=', G, '\n')
for (r in 1:G) print(round(apply(theta.pomtc.true[r,,,], c(1,2), mean),2))

# Sampling y
for (j in 2:(p+1)) for (i in 1:n) 
  y[i,j]=sample(1:q,1,prob=theta.pomtc.true[r.true[i],y[i,j-1],,(j-1)],replace=T)
write.csv(y,"y.sim.csv", row.names=F) 

# First observation is not observable (cluster membership is random at start)
y.mat=y[,-1]

#############################################################################
######################       (2) Run model            #######################
# Read data
cat(" ####################################################################",'\n')
cat('DATASET:', thisdata, '\n')
str(y.mat)
y.mat=matrix(unlist(y.mat),nrow=nrow(y.mat),byrow=F)
mydata = list(y=y.mat)
n=dim(mydata$y)[1]
p=dim(mydata$y)[2]
q=length(unique(as.numeric(mydata$y)))

# MLE Regression for inits
pom.mle <- polr(as.factor(y) ~ 1, data=mydata,method='logistic')
summary(pom.mle)
mu.mle=pom.mle$zeta

############################    PROPOSALS    ##############################
#alpha
proposal.alpha=sqrt(0.1^2) 
proposal.var.alpha=sqrt(log(1.01))  
#beta 
proposal.beta=sqrt(0.1^2)
proposal.var.beta=sqrt(log(1.01)) 
#gamma
proposal.gamma=sqrt(0.1^2)
proposal.var.gamma=sqrt(log(1.01))  
# mu
proposal.mu=0.1  # 0.075
proposal.var.mu=sqrt(log(1.01))  
# pi
proposal.sd.pi=0.1

# Parameters for proposals
qpar <- list(proposal.mu=proposal.mu, proposal.var.mu=proposal.var.mu,
             proposal.alpha=proposal.alpha, proposal.var.alpha=proposal.var.alpha,
             proposal.beta=proposal.beta, proposal.var.beta=proposal.var.beta,
             proposal.gamma=proposal.gamma, proposal.var.gamma=proposal.var.gamma,
             proposal.sd.pi=proposal.sd.pi)
cat("My proposals", '\n')
print(qpar)
cat("DATA:", thisdata, " n=",n,' p=',p," q=",q, sep='', '\n')

# Printing out true values
param.true = list(mu=mu.true, sigma2.mu=3,  # this depends on the reparametrization chosen for mu!
                  alpha=alpha.true, sigma2.alpha=sigma2.gamma.true,
                  beta=beta.true, sigma2.beta=sigma2.beta.true,
                  gamma=gamma.true, sigma2.gamma=sigma2.gamma.true,
                  pi=pi.true )
state.true <- init.state(param.true, mydata, fpar, qpar)
cat("TRUE mu=",round(mu.true,2)," alpha=",round(alpha.true,2)," sigma2.beta=",round(sigma2.beta.true,2), " pi=",round(pi.true,2), " log.like.true=",  round(state.true$log.like,1)," log.like.prior=",  round(state.true$log.prior,1), sep=' ','\n')
cat('gamma=',round(gamma.true,2), '\n')
cat('beta')
print(round(beta.true,2))

# Parameters for Priors
vars=c("mu", "sigma2.mu","alpha", "sigma2.alpha", "beta", "sigma2.beta","gamma","sigma2.gamma", "pi")
# beta_rk'
if (R==1) npars=(q-1)+1+ (q-1)+1+(p-2)+1 else 
  npars=(q-1)+1+(R-1)+1+ (R)*(q-1)+1+(R-1)+(p-2)+1 # number of independent parameters
nfollow <-  (q-1)+1+R+1+R*(q+1)+R+p+1+2+9             # all pars + like, post + AR
# saving all to a list
fpar <- list(a=mya, b=myb, phi=3/2*rep(1,R), n=n, p=p ,q=q, R=R, npars=npars, nfollow=nfollow, vars=vars)
# printing out 
cat('Parameters for IG prior: shape=', fpar$a, ' rate=', fpar$b,' mode=',fpar$b/(fpar$a+1), 'mean', round(fpar$b/(fpar$a-1),1),'\n')

# Starting values
Inits=list()
POM.rcc.bayes =list()
# Starting values for pi and alpha based on k-means
kmeans.data=kmeans(y.mat,R,nstart=50)
pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
alpha.kmeans=apply(kmeans.data$centers,1,mean)
set.seed(myseed)

# beta_rk'=0, gamma_j=0
Inits[[1]] = list(mu=mu.mle-mu.mle[1], sigma2.mu=1,
                  alpha=alpha.kmeans, sigma2.alpha=1,
                  #beta=rep(0,q), sigma2.beta=1,
                  beta=matrix(rep(0,R*q),nrow=R), sigma2.beta=rep(0.25,R),
                  gamma=rep(0,p), sigma2.gamma=1,
                  pi=pi.kmeans)
# alpha=0
Inits[[2]] = list(mu=mu.mle-mu.mle[1], sigma2.mu=1,
                  alpha=rep(0,R), sigma2.alpha=1,
                  beta=matrix(rep(0,R*q),nrow=R), sigma2.beta=rep(0.25,R),
                  gamma=rep(0,p), sigma2.gamma=1,
                  pi=pi.kmeans)
# all above plus same pi
Inits[[3]] = list(mu=mu.mle-mu.mle[1], sigma2.mu=1,
                  alpha=rep(0,R), sigma2.alpha=1,
                  beta=matrix(rep(0,R*q),nrow=R), sigma2.beta=rep(0.25,R),
                  gamma=rep(0,p), sigma2.gamma=1,
                  pi=rep(1/R,R))

# Runing several chains using lappy
t0 <- proc.time()[3]
ltemp=lapply(1:ninits,pomtc,mydata, fpar, qpar)
ta <- round((proc.time()[3] - t0)/60,2)

cat("\n \n pomtc run time is", ta, "mins for lapply with",ninits," inits for n,p,q,R",n,p,q,R," \n")

# Saving results
cat("  Saving R session with MCMC results \n")
name=paste("pomtc_",thisdata,'_n',n,"_p",p,"_q",q,"_R",R ,"_ninits",ninits,"_t",ta, "min.RData", sep='')
save.image(name)


#############################################################################
#####################     (3) Output postprocessing   #######################
#####################     relabel MCMC chains, etc    #######################
# Post processing here
cat(" \n ####################################################   ")
cat("  Starting Post processing MCMC chains here \n \n")

cat("TRUE mu=",round(mu.true,2)," alpha=",round(alpha.true,2)," sigma2.beta=",round(sigma2.beta.true,2), " pi=",round(pi.true,2), " log.like.true=",  round(state.true$log.like,1)," log.like.prior=",  round(state.true$log.prior,1), sep=' ','\n')
cat('gamma=',round(gamma.true,2), '\n')
cat('beta')
print(round(beta.true,2))

# Converting chains as mcmc objects
stopifnot(ninits<=4)
for (i in 1:ninits) ltemp[[i]]=mcmc(ltemp[[i]])
allchains=mcmc.list(ltemp[[1]])
if (ninits==2) allchains=mcmc.list(ltemp[[1]],ltemp[[2]])
if (ninits==3) allchains=mcmc.list(ltemp[[1]],ltemp[[2]],ltemp[[3]])
if (ninits==4) allchains=mcmc.list(ltemp[[1]],ltemp[[2]],ltemp[[3]],ltemp[[4]])

# plot for all chains
name=paste('pomtc_',thisdata,'_n',n,"_p",p, "_q",q,"_R",R,'_ninits',ninits,".pdf",sep="")
pdf(name,width=6, height=8,pointsize=10)
plot(allchains)
dev.off()

cat(" \n ####################################################   ")
cat(" Relabelling chains using Stephen's algorithm \n")
chain=ltemp[[1]]
if (ninits>1) for (init in 2:ninits) chain=rbind(chain, ltemp[[1]] )
D=dim(chain)[1] 

betas=paste("beta",1,'.',1:(q),sep='')
if (R>1) for (r in 2:R) betas=c(betas,paste("beta",r,'.',1:(q),sep=''))
deviance=-2*chain[,"llike"]
log.like=chain[,"llike"]
d.bar = mean(deviance)
pv=0.5*var(deviance)
dic.pv=d.bar+pv
# plug-in values
mu.plug = apply(chain[ ,paste("mu",1:(fpar$q-1),sep='')],2,mean)
alpha.plug = apply(matrix(chain[ ,paste("alpha",1:(fpar$R),sep='')],ncol=fpar$R),2,mean) # matrix to make APPLY
beta.plug = matrix(apply(chain[ , betas],2,mean), byrow = T, nrow=fpar$R)
gamma.plug = apply(chain[ ,paste("gamma",1:(fpar$p),sep='')],2,mean)
pi.plug= apply(matrix(chain[, paste("pi",1:(fpar$R),sep='')],ncol=fpar$R),2,mean) # matrix to make APPLY
theta.plug = theta.pomtc(mu.plug, alpha.plug, beta.plug, gamma.plug, fpar)
Z.plug=Z.rc.tt(theta.plug, pi.plug,fpar,mydata)
log.like.plug=sum(log(apply(exp(Z.plug),1,sum)))
cat("Plug-in (means): \n")
cat("mu=",round(mu.plug,2)," alpha=",round(alpha.plug,2), " pi=",round(pi.plug,2), 'log.like=', round(log.like.plug,2), '\n')
print(round(gamma.plug,2))
print(round(beta.plug,2))
d.plug = -2*log.like.plug
p.dic = d.bar - d.plug
dic = d.bar + p.dic
cat("d.bar=",round(d.bar,2), " p=",round(p.dic,2), " dic=", round(dic,2), " pv=",round(pv,2), " dic.pv=", round(dic.pv,2),'\n')

# WAIC
Z.temp=matrix(99,n,fpar$R) 
Z.mcmc=array(NA,dim=c(n,D))
cat('Waic calculations','\n')
for (d in 1:D) {
  mu = chain[d,paste("mu",1:(fpar$q-1),sep='')]
  alpha = chain[d ,paste("alpha",1:(fpar$R),sep='')]
  beta = matrix(chain[d, betas], byrow=T, nrow=fpar$R)
  gamma = chain[d,paste("gamma",1:(fpar$p),sep='')]
  pi= chain[d, paste("pi",1:(fpar$R),sep='')]
  theta = theta.pomtc(mu,alpha,beta,gamma,fpar)
  Z = exp(Z.rc.tt(theta,pi,fpar,mydata)) 
  Z.mcmc[,d]= log(apply(Z,1,sum))
  if( d/100==floor(d/100)) cat('*')
}
cat('\n')
v=sum(apply(Z.mcmc,1,var))
waic.g=d.bar+2*v

# Z is in logs
p.waic1 = 2* sum( log(apply(exp(Z.mcmc),1,mean)) - apply(Z.mcmc,1,mean) )
p.waic2= sum(apply(Z.mcmc,1,var))
lppd=-2*sum(log(apply(exp(Z.mcmc),1,mean)))
waic1=lppd+2*p.waic1
waic2=lppd+2*p.waic2
cat("lppd=",round(lppd,1)," p.waic1=",round(p.waic1,1)," waic1=", round(waic1,1), " p.waic2=",round(p.waic2,1)," waic2=", round(waic2,1),'\n')
if (R==1) {
  npars.like= (q-1)+ (q-1) # no alpha, no pi
} else {
  npars.like= (q-1)+ (R-1)+ R*(q-1)+ R
}

cat("###################    Total Draws=", D," from ninits=", ninits, "each" , dim(ltemp[[1]])[1], "  ###################", '\n')

if (relabel==1 & fpar$R>1) {
  cat("RELABELLING: relabel=", relabel, " R=", fpar$R, '\n')
  #########################################################################
  # Label Switching (Stephen's algorithm)
  cat("Relabelling chains", '\n')
  chain1=chain # chain to store temp values
  chain2=chain # new chain to store final relabeled values
  iter=1
  nperm= factorial(fpar$R)
  no.convergence=TRUE
  while ( (iter==1) | no.convergence  ) {
    cat("Iter=", iter, '\n')
    chain1 <- chain2
    ###################################################
    ###############                 Step1
    # Calculate Q=E[Z.mcmc], ie qir is mean of Zir_mcmc (over mcmc draws)
    cat(' calculating Q \n')
    Z.mcmc=array(NA,dim=c(fpar$n,fpar$R, D))
    mu = chain1[,paste("mu",1:(fpar$q-1),sep='')]
    alpha = chain1[ ,paste("alpha",1:(fpar$R),sep='')]
    beta = chain1[ ,betas]
    gamma = chain1[ ,paste("gamma",1:(fpar$p),sep='')]
    pi= chain1[, paste("pi",1:(fpar$R),sep='')]
    sigma2.beta = chain1[ ,paste("sigma2.beta",1:(fpar$R),sep='')]
    # Objects to save coclustering matrix
    r.mcmc=matrix(NA,n,D)
    entropy.mcmc=matrix(NA,D)
    coclustering.mat.mean=matrix(0,n,n)
    Z.mcmc.mean=matrix(0,n,R)
    for (d in 1:D) {
      theta = theta.pomtc(mu[d,],alpha[d,], beta[d,], gamma[d,],fpar)
      Z =Z.rc.tt(theta,pi[d,],fpar,mydata) # already in logs, only num of Z (not yet a posterior prob)
      Z=exp(Z)/apply(exp(Z),1,sum)
      Z.mcmc[,,d]=Z
      Z.mcmc.mean= Z + Z.mcmc.mean
      r.mcmc[,d]=group(Z)
      entropy.mcmc[d]=entropy(Z)
      temp.mat=matrix(as.numeric(r.mcmc[,d] %*% t(r.mcmc[,d]) %in% (1:q)^2),n)
      coclustering.mat.mean=coclustering.mat.mean+temp.mat
    }
    Q=apply(Z.mcmc,c(1,2),mean)
    ###################################################
    ###############                 Step2
    cat(' Permutation with lowest loss')
    Loss=array(-99,dim=c(D, nperm))
    for (d in 1:D) {
      # Permuting the indexes not values
      allperms = permutations(n=fpar$R , r=fpar$R, v=1:fpar$R)
      alpha.per=allperms
      pi.per=allperms
      sigma2.beta.per=allperms
      beta.per=array(-99, dim=c(nperm, fpar$R ,q)) ## beta_rk' is a R(p-1) matrix
      beta.d = t( matrix(beta[d,], ncol=fpar$R) )
      for (per in 1:nperm) {
        for (g in 1:fpar$R) {
          alpha.per[per,g]= alpha[d,][allperms[per,g]]
          pi.per[per,g]= pi[d,][allperms[per,g]]
          sigma2.beta.per[per,g]= sigma2.beta[d,][allperms[per,g]]
          beta.per[per, , ][g,] = beta.d[allperms[per,g],]
        }
      }
      mu.per=t(matrix(rep(mu[d,], nperm), ncol=nperm))
      # Creating an array to hold theta for each permutation
      theta.per=array(99, dim=c(nperm,fpar$R,fpar$q,fpar$q,fpar$p))
      Z.per=log.Z.per=array(99, dim=c(nperm,fpar$n,fpar$R))
      # No label switching in gamma
      for (per in 1:nperm) {
        theta.per[per,,,,] = theta.pomtc(mu.per[per,], alpha.per[per,], beta[d,], gamma[d,], fpar)
        log.Z.per[per,,] = Z.rc.tt(theta.per[per,,,,], pi.per[per,], fpar, mydata)
        Z.per[per,,] = exp(log.Z.per[per,,])/apply(exp(log.Z.per[per,,]),1,sum)
        Loss[d,per] = sum(Z.per[per,,]*(log.Z.per[per,,]-log(Q)))
      }
      minper = which.min(Loss[d,])  #maxperm = which.max(Loss[i,])
      # mu and gamma do not change
      chain2[d, paste("mu",1:(fpar$q-1),sep='') ] = chain1[d,paste("mu",1:(fpar$q-1),sep='')]
      chain2[d, paste("gamma",1:(fpar$p),sep='') ]  = chain1[d ,paste("gamma",1:(fpar$p),sep='')]
      # changing parameters affected by label switching
      chain2[d , paste("alpha",1:(fpar$R),sep='') ]=  alpha.per[minper,]
      chain2[d , betas ] = as.numeric(t(beta.per[minper, ,]))
      chain2[d, paste("sigma2.beta",1:(fpar$R),sep='')] =  sigma2.beta.per[minper,]
      chain2[d, paste("pi",1:(fpar$R),sep='')] =  pi.per[minper,]
      if( d/100==floor(d/100) ) cat('*')
    }
    iter=iter+1
    no.convergence <- any(abs(apply(chain1[ ,paste("alpha",1:(fpar$R),sep='')],2,mean)-apply(chain2[ ,paste("alpha",1:(fpar$R),sep='')],2,mean))> tol)
    no.convergence
    cat(" chain 1 mean alpha=",round(apply(chain1[ ,paste("alpha",1:(fpar$R),sep='')],2,mean),3),"\n",
        "chain 2 mean alpha=",round(apply(chain2[ ,paste("alpha",1:(fpar$R),sep='')],2,mean),3),"\n",
        'Convergence:', ! no.convergence,'Differences alpha: ', round(abs(apply(chain1[ ,paste("alpha",1:(fpar$R),sep='')],2,mean)-apply(chain2[ ,paste("alpha",1:(fpar$R),sep='')],2,mean)),3),'\n')
  }
  #end Stephen's algorithm for relabelling
  
  #### Graphs and GOF  with relabelled chain2  ####
  chain1 = chain # saving original relabelled chain in chain1
  chain = chain2
  ##########################
  cat(" chain mean alpha=",round(apply(chain[ ,paste("alpha",1:(fpar$R),sep='')],2,mean),3),"\n")
  
  # Graphs of relabelled chains
  name=paste('pomtc_',thisdata,'_n',n,"_p",p, "_q",q,"_R",R,'allrelabelled_ninits',ninits,".pdf",sep="")
  pdf(name,width=6, height=8,pointsize=10)
  plot( mcmc.list( mcmc(chain) ),ask=F)
  dev.off()
  #########
  # DIC
  deviance=-2*chain[,"llike"]
  log.like=chain[,"llike"]
  d.bar = mean(deviance)
  pv=0.5*var(deviance)
  dic.pv=d.bar+pv
  
  # Posterior median
  cat("Plug-in is median")
  mu.plug = apply(chain[ ,paste("mu",1:(fpar$q-1),sep='')],2,median)
  alpha.plug = apply(matrix(chain[ ,paste("alpha",1:(fpar$R),sep='')],ncol=fpar$R),2,median)
  beta.plug = matrix(apply(chain[ , betas],2,median), byrow = T, nrow=fpar$R)
  gamma.plug = apply(chain[ ,paste("gamma",1:(fpar$p),sep='')],2,median)
  pi.plug= apply(matrix(chain[, paste("pi",1:(fpar$R),sep='')],ncol=fpar$R),2,median)
  
  #theta,Z and like
  theta.plug = theta.pomtc(mu.plug, alpha.plug, beta.plug, gamma.plug, fpar)
  Z.plug=Z.rc.tt(theta.plug, pi.plug,fpar,mydata)
  log.like.plug=sum(log(apply(exp(Z.plug),1,sum)))
  cat("Plug-in (means): \n")
  cat("mu=",round(mu.plug,2)," alpha=",round(alpha.plug,2), " pi=",round(pi.plug,2), 'log.like=', round(log.like.plug,2), '\n')
  print(round(gamma.plug,2))
  print(round(beta.plug,2))
  d.plug = -2*log.like.plug
  p.dic = d.bar - d.plug
  dic = d.bar + p.dic
  cat("d.bar=",round(d.bar,2), " p=",round(p.dic,2), " dic=", round(dic,2), " pv=",round(pv,2), " dic.pv=", round(dic.pv,2), '\n')
  ###########
}

# Re-ordering clusters by increasing alpha
# true values
myorder = order(alpha.true)
cat( "Order in true values \n", myorder)
if ( any(myorder!=1:R) ) {
  cat( "Re-ordering true values \n", myorder)
  alpha.true = alpha.true[myorder]
  pi.true = pi.true[myorder]
  sigma2.beta.true = sigma2.beta.true[myorder]
  temp=beta.true
  for (r in 1:R) beta.true[r,]=temp[myorder[r],] 
} else cat("No need to re-order true values by incresing alpha \n")
# saving all true values 
truesim = c(mu.true, NA,
            alpha.true, NA,
            as.numeric(t(beta.true)), sigma2.beta.true,
            gamma.true,sigma2.gamma.true,
            pi.true,
            state.true$log.like,
            state.true$log.like+state.true$log.prior)
# raw MCMC results
convergence.results = cbind(mean=apply(chain,2,mean),se=apply(chain,2,sd), HPDinterval(mcmc(chain)))

# output a table to save all ordered results
output = truesim
names(output)=colnames(chain)
alpha.post = apply(chain[,paste("alpha",1:R,sep='')],2, mean) 
myorder=order(alpha.post)
cat( "Order in posterior values \n", myorder, "\n")

for (stat in colnames(convergence.results)){
  cat(" Re-ordering posterior parameters for ", stat, " mixture labels", myorder,'\n')
  
  if ( any(myorder!=1:R) ) {
    # Getting values from chain
    post.stat = convergence.results[,stat]
    # each parameter
    mu.post =  post.stat[paste("mu",1:(q-1),sep='')]
    sigma2.mu.post = post.stat["sigma2.mu"]
    alpha.post = post.stat[paste("alpha",1:R,sep='')]
    sigma2.alpha.post = post.stat["sigma2.alpha"]
    pi.post = post.stat[paste("pi",1:R,sep='')]
    sigma2.beta.post = post.stat[paste("sigma2.beta",1:R,sep='')]
    betas = paste("beta",1,".",1:(q),sep='')
    for (r in 2:R) betas=c(betas,paste("beta",r,".",1:(q),sep=''))
    beta.post = matrix( post.stat[betas], byrow = T, nrow = R)
    gamma.post = post.stat[paste("gamma",1:p,sep='')]
    sigma2.gamma.post = post.stat["sigma2.gamma"]
    llike = post.stat["llike"]
    lpost = post.stat["lpost"]
    # Orderning
    alpha.post = alpha.post[myorder]
    pi.post = pi.post[myorder]
    sigma2.beta.post = sigma2.beta.post[myorder]
    temp=beta.post
    for (r in 1:R) beta.post[r,]=temp[myorder[r],] 
    
    # re-ordered parameters
    post.o = c(mu.post, NA,
               alpha.post, NA,
               as.numeric(t(beta.post)), sigma2.beta.post,
               gamma.post,sigma2.gamma.post,
               pi.post,
               llike,
               lpost)
    output = cbind(output, post.o)
  } else cat("No need to order posterior values by increasing alpha  \n")
}  
colnames(output) = c("true", colnames(convergence.results))
print(round(output,2))
write.csv(round(output,2), paste("sim.resuts.n",n,".p",p,".q",q,".G",G,".csv",sep=""))

#   Latex format for convergence results
library(xtable)
allnames <- c(paste("$\\mu_",1:(q-1),"$",sep=''),"$\\sigma^2_{\\mu}$")
allnames <- c(allnames, paste("$\\alpha_",1:(R),"$",sep=''), "$\\sigma^2_{\\alpha}$")
betas_latex=paste("$\\beta_{",1,1:(q),"}$",sep='')
betas_latex
if (R>1) for (r in 2:R) betas_latex=c(betas_latex,paste("$\\beta_{",r,1:(q),"}$",sep=''))
betas_latex
allnames <- c(allnames, betas_latex)
allnames <- c(allnames, paste("$\\sigma^2_{\\beta",1:R,"}$",sep=""))
allnames <- c(allnames, paste("$\\gamma_{",1:p,"}$",sep=''), "$\\sigma^2_{\\gamma}$")
allnames <- c(allnames, paste("$\\pi_",1:(R),"$",sep=''))
allnames <- c(allnames, "log-like","log-post")
rownames(output)=allnames
output.table = xtable(output, caption=paste('pomtt Convergence results for',thisdata, 'data'))
digits(output.table) = 2
print.xtable(output.table, sanitize.text.function = function(x) x)

# Saving results with relabelled chains
cat("  Saving R session with results again \n")
name=paste("pomtc_",thisdata,'_n',n,"_p",p,"_q",q,"_G",G ,"_ninits",ninits,"_t",ta, "min.RData", sep='')
save.image(name)
cat(" \n #######################        DONE!!!!     #######################    \n ")