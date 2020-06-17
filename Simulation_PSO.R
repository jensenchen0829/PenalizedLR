library(pso)
library(GA)
n<-1000
mo<-4 #the prior density of beta is Normal
p<-19
simUnmatched = function(n,p,scale=FALSE){
  # n is total sample size, beta1 is value of parameter of interest,
  # p is number of nuisance covariates.wqs
  ConCaseRatio = 4 # assuming 4:1 con:case ratio
  ncase = n/(ConCaseRatio+1); ncon=ncase*ConCaseRatio
  beta = log(rf(p+1,mo/2,mo/2)) # p nuisance params of value 1
  ncov = p+1
  # Simulate cases and controls
  conX = caseX = NULL
  for(i in 1:ncov) {
    conX = cbind(conX,rnorm(ncon,mean=0,sd=1))
    caseX = cbind(caseX,rnorm(ncase,mean=beta[i],sd=1))
  }
  X = rbind(caseX,conX)
  colnames(X) = paste0("x",1:ncov); rownames(X) = NULL
  case = c(rep(1,ncase),rep(0,ncon))
  return(data.frame(case,X))
}


psoLA<-function(alpha_star.V,m,n.rounds){
  tracer<-matrix(0,nrow=1,ncol=p+4)
  ftracer=0
  i=1
  for (i in 1:n_rounds){
    beta_max<-numeric(p+1)
    for (di in 1:(p+1)){
      X<-XM[,di]
      alpha_star<-alpha_star.V[di]
      dlogPenalisedL<-function(beta){
        sum(X*y-(X*exp(alpha_star+beta*as.numeric(X))/(1+exp(alpha_star+beta*as.numeric(X)))))-m/2+m*exp(-beta)/(1+exp(-beta))
      }
      beta_max[di]<-uniroot(dlogPenalisedL, c(-20,20))$root
    }
    multi.dimen.logLP_betamax<-function(alpha_star0m0){
      m0=alpha_star0m0[p+2]
      ll<-0
      for (di in 1:(p+1)){
        alpha_star0<-alpha_star0m0[di]
        X<-XM[,di]
        temp1<-sum(X^2*exp(alpha_star0+beta_max[di]*X)/(1+exp(alpha_star0+beta_max[di]*X)))-sum(X^2*(exp(alpha_star0+beta_max[di]*X)/(1+exp(alpha_star0+beta_max[di]*X)))^2)
        temp2<-exp(-beta_max[di])/(1+exp(-beta_max[di]))-(exp(-beta_max[di])/(1+exp(-beta_max[di])))^2
        c=temp1+temp2
        LP_di<-sum(y*(alpha_star0+beta_max[di]*as.numeric(X))-log(1+exp(alpha_star0+beta_max[di]*as.numeric(X))))-log(beta(m0/2,m0/2))-m0/2*beta_max[di]-m0*log(1+exp(-beta_max[di]))+0.5*log(c)
        ll<-ll+LP_di
      }
      ll
    }
    pso.result<-psoptim(par=c(alpha_star.V,m),fn=multi.dimen.logLP_betamax,lower=c(rep(-20,p+1), 0), upper = c(rep(10,p+1),20),control=list(trace=100,fnscale=-1,maxit=2000,maxit.stagnate=30,s=45,type="SPSO2011"))
    if(abs(pso.result$value-ftracer)>=0.001*abs(ftracer)){
      alpha_star.V=pso.result$par[1:(p+1)]
      m=pso.result$par[p+2]
      tracer<-rbind(tracer,c(as.integer(i),alpha_star.V,m,pso.result$value))
      ftracer<-pso.result$value
      print(tracer[i+1,])
    }else{
      break
    }
  }
  tracer
}


###remember to change s

  simdata<-simUnmatched(n,p)
  y<-simdata$case
  XM<-as.matrix(simdata[,-1])
  #glm1<-glm(y~XM, family=binomial)
  ##Initialization
  ini.alpha<-rep(unname(-3),p+1)
  ini.m<-6
  n_rounds<-400
  tracer1<-psoLA(ini.alpha,ini.m)

