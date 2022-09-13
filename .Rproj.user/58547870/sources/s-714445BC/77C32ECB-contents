#' @rdname robustbetareg
#'
#' @param x,z  numeric regressor matrix for mean and precision model respectively, defaulting to an intercept only.
#'
# Robust Estimation - LSMLE
LSMLE.fit=function(y,x,z,alpha=NULL,link="logit",link.phi="log",control=robustbetareg.control(...),...)
{
  #options(warn = 2) #Convert warnings in errors
  result=theta=list()
  #Arguments Checking
  alpha.optimal=control$alpha.optimal
  start_theta=control$start
  if(!is.null(alpha)){alpha.optimal=FALSE}
  if(is.null(alpha) & alpha.optimal==FALSE){alpha=0}
  linkobj=set.link(link.mu = link,link.phi = link.phi)
  k=ncol(x)
  m=ncol(z)
  if(alpha.optimal)
  {
    return(Opt.Tuning.LSMLE(y,x,z,link,link.phi,control))
  }
  if(is.null(start_theta))
  {
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
    theta$x=start_theta
    theta$converged=mle$converged
  }
  #Point Estimation
  #browser()
  q=1-alpha
  check=TRUE
  theta=tryCatch(optim(par=start_theta,fn=L_alpha_R,gr=Psi_LSMLE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1,maxit=10000)),error=function(e){
    theta$msg<-e$message;check<<-F;return(theta)})
  if(check){
    if(theta$convergence==0){
      theta$converged=T
      theta$x=theta$par
    }else{
      theta$converged=F
      theta$x=start_theta
    }
  }else{
    theta$converged=F
    theta$x=start_theta
  }

  #Predict values
  beta=theta$x[1:k]
  gamma=theta$x[seq.int(length.out = m) + k]
  coefficients = list(mean = beta, precision = gamma)
  eta=x%*%beta
  mu_hat=linkobj$linkfun.mu$inv.link(eta)
  phi_hat=linkobj$linkfun.phi$inv.link(z%*%gamma)

  #Expected Standard Error
  M.LSMLE=tryCatch(LSMLE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj),error=function(e) NULL)
  if(is.null(M.LSMLE))
  {
    vcov=std.error.LSMLE=NaN
  }else{
    vcov=M.LSMLE$Cov
    std.error.LSMLE=M.LSMLE$Std.Error
  }

  #Register of output values
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients#Coefficients Regression
  result$vcov=vcov#Expected Covariance Matrix
  pseudor2 = if (var(eta) * var(qlogis(y)) <= 0){NA}
  else{cor(eta, linkobj$linkfun.mu$linkfun(y))^2}
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.LSMLE)))#Convergence
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)#Fitted Values
  result$start=start_theta #Started Point
  result$weights=(degbeta(y_star,mu_hat,phi_hat/q))^(alpha)#Weights
  result$Tuning=alpha#Tuning
  result$residuals=sweighted2_res(mu_hat,phi_hat,y=y,X=x,linkobj = linkobj)#Standard Residual Report
  result$n=length(mu_hat)#Sample Size
  result$link=link#Mean Link Function
  result$link.phi=link.phi#Precision Link Function
  result$Optimal.Tuning=alpha.optimal#Flag - Data-driven algorithm tuning selector
  result$pseudo.r.squared=pseudor2
  result$control=control#Extra options package
  result$call=match.call()
  result$method="LSMLE"
  if(any(is.na(std.error.LSMLE)))
  {
    str1="Standard-Error is unvailable"
  }
  if(!is.null(theta$msg))
  {
    str2=theta$msg
  }
  if(!is.null(str1)||!is.null(str2))
  {
    result$Message=c(str1,str2)
  }
  if(!any(is.na(std.error.LSMLE)))
  {#Standard Error
    names(std.error.LSMLE)<-c(colnames(x),colnames(z))
    se.beta=std.error.LSMLE[1:k]
    se.gamma=std.error.LSMLE[1:m+k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error=coef.std.error
  }
  class(result)="robustbetareg"
  return(result)
}


###################################################################################################################


# Auto Selecting tuning parameter algorithm
Opt.Tuning.LSMLE=function(y,x,z,link,link.phi,control)
{
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LSMLE.list=LSMLE.par=list()
  zq.t=NULL
  alpha_tuning=seq(0,0.6,0.02)
  K=length(alpha_tuning)
  M1=11
  M=control$M
  L=control$L
  n=length(y)
  unstable=F
  sqv.unstable=T

  est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
  if(is.null(est.log.lik))
  {
    est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi, control=betareg.control(start=Initial.points(y,x,z)))),error=function(e) NULL)
  }
  if(!is.null(est.log.lik)){
    Est.param=do.call("c",est.log.lik$coefficients)
    names(Est.param)=c(colnames(x),colnames(z))
  }else{
    Est.param=Initial.points(y,x,z)
  }
  if(!is.null(control$start)){
    Est.param=control$start
  }
  ponto.inicial.robst=Initial.points(y,x,z)
  ponto.inicial.temp=Est.param
  names(ponto.inicial.robst)=c(colnames(x),colnames(z))
  p=length(Est.param)
  for(k in 1:M1)#First M1 attempts of tuning parameters
  {
    control$start=ponto.inicial.temp
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=Est.param
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=c(LSMLE.par$coefficients$mean,LSMLE.par$coefficients$precision)
    }
    if(!LSMLE.par$converged || is.null(LSMLE.par) || any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LSMLE.list[[k]]<-LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
  }
  #sqv=as.numeric(SQV_Cpp(zq.t,n,p))
  sqv=as.numeric(SQV(zq.t,n,p))
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    LSMLE.par.star<-LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  if(unstable)#Lack of stability
  {
    LSMLE.par.star=LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    LSMLE.par.star$message="Lack of stability"
    return(LSMLE.par.star)
  }
  if(alpha.ind<8){#Which Tuning satisfy the condition os stability
    LSMLE.par.star<-LSMLE.list[[alpha.ind+1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }

  reached=FALSE
  k=M1+1
  while(sqv.unstable & !reached)#Seek within the next grid of tuning
  {
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=c(LSMLE.par$coefficients$mean,LSMLE.par$coefficients$precision)
    }
    if(any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      unstable=T
      break
    }
    LSMLE.list[[k]]=LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))

    sqv=as.numeric(SQV(zq.t,n,p))
    sqv.test=sqv[(k-M):(k-1)]
    if(all(sqv.test<=L))
    {
      sqv.unstable=F
    }
    k=k+1
    if(k>=K)#Condition of Step 6
    {
      reached=TRUE
    }
  }
  if(reached)
  {
    k=suppressWarnings(max(1,min(which(rollapply(sqv<L,M,sum)==M)))+M+1)
  }
  if(k>=K || unstable)
  {
    LSMLE.par.star=LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    LSMLE.par.star$message="Lack of stability"
  }else{
    LSMLE.par.star=LSMLE.list[[(k-1-M)]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
  }
  return(LSMLE.par.star)
}


#Objective Function - LSMLE Log-likelihood
L_alpha_R=function(theta,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]
  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)
  phi_q = phi_hat/q
  y_star=log(y)-log(1-y)
  f_q_star = degbeta(y_star,mu_hat,phi_q)
  if(alpha==0){
    L_q=sum(log(f_q_star))
  }else{
    L_q=sum((f_q_star^(alpha)-1)/alpha)
  }
  return(L_q)
}

#Gradiente: Psi_beta LSMLE
Psi_Beta_LSMLE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  phi_q=phi_hat/q

  aq=mu_hat*phi_q
  bq=(1-mu_hat)*phi_q

  y_star=log(y)-log(1-y)
  mu_star=digamma(aq)-digamma(bq)
  Tb = (link.model$link.mu$d.linkfun(mu_hat))^(-1)
  f_q_star = (degbeta(y_star,mu_hat,phi_q))^(alpha)

  result=t(X)%*%diag(phi_q*Tb*f_q_star)%*%(y_star-mu_star)
  return(result)
}

#Gradiente: Psi_gamma LSMLE
Psi_Gamma_LSMLE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  phi_q = phi_hat/q
  aq = mu_hat * phi_q
  bq = (1-mu_hat)*phi_q
  y_dagger = log(1-y)
  y_star = log(y)-y_dagger
  mu_star = digamma(aq)-digamma(bq)
  mu_dagger = digamma(bq)-digamma(phi_q)
  Tg = (link.model$linkfun.phi$d.linkfun(phi_hat))^1

  f_q_star = (degbeta(y_star,mu_hat,phi_q))^(alpha)
  eta = mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger)
  result=t(Z)%*%diag(Tg%*%f_q_star/q)%*%eta
  return(result)
}

#Gradiente: Psi LSMLE
Psi_LSMLE=function(theta,y,X,Z,alpha,link_mu,link_phi){
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]
  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  psi_beta = Psi_Beta_LSMLE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)
  psi_gamma = Psi_Gamma_LSMLE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)

  return(c(psi_beta,psi_gamma))
}

#Sandwich Matrix - LSMLE
LSMLE_Cov_Matrix=function(mu,phi,X,Z,alpha,linkobj)
{
  n=length(mu)
  q=1-alpha
  #
  a.0=mu*phi
  b.0=(1-mu)*phi
  #
  phi.k0=phi/(1-alpha)
  a.k0=mu*phi.k0
  b.k0=(1-mu)*phi.k0
  #
  phi.k1=phi*(1+alpha)/(1-alpha)
  a.k1=mu*phi.k1
  b.k1=(1-mu)*phi.k1

  C.k0=diag(exp(q*lbeta(a.k0,b.k0)-lbeta(a.0,b.0)))
  C.k1=diag(exp(lbeta(a.k1,b.k1)-lbeta(a.0,b.0)-2*alpha*lbeta(a.k0,b.k0)))

  mu_star.k0=digamma(a.k0)-digamma(b.k0)
  mu_dagger.k0=digamma(b.k0)-digamma(phi.k0)
  mu_star.k1=digamma(a.k1)-digamma(b.k1)
  mu_dagger.k1=digamma(b.k1)-digamma(phi.k1)

  #mean link function
  Tb=diag((linkobj$linkfun.mu$d.linkfun(mu))^(-1))
  #precision link funtion
  Tg=diag((linkobj$linkfun.phi$d.linkfun(phi))^(-1))

  Phi=diag(phi)
  Q.inv=diag(n)/q

  Lambda_mu_mu=diag(trigamma(a.k0)+trigamma(b.k0))
  Lambda_mu_phi=diag(mu*trigamma(a.k0)-(1-mu)*trigamma(b.k0))
  Lambda_phi_phi=diag(mu^2*trigamma(a.k0)+(1-mu)^2*trigamma(b.k0)-trigamma(phi.k0))

  Lambda_beta_beta=t(X)%*%Tb%*%Q.inv%*%(Phi^2)%*%C.k0%*%Lambda_mu_mu%*%Tb%*%X
  Lambda_beta_gamma=t(X)%*%Tb%*%Q.inv%*%Phi%*%C.k0%*%Lambda_mu_phi%*%Tg%*%Z
  Lambda_gamma_gamma=t(Z)%*%Tg%*%Q.inv%*%C.k0%*%Lambda_phi_phi%*%Tg%*%Z

  Lambda=-1*rbind(cbind(Lambda_beta_beta,Lambda_beta_gamma),cbind(t(Lambda_beta_gamma),Lambda_gamma_gamma))

  Sigma_mu_mu=diag(trigamma(a.k1)+trigamma(b.k1)+(mu_star.k1-mu_star.k0)^2)
  Sigma_mu_phi=diag(mu*trigamma(a.k1)-(1-mu)*trigamma(b.k1)+mu*(mu_star.k1-mu_star.k0)^2+(mu_star.k1-mu_star.k0)*(mu_dagger.k1-mu_dagger.k0))
  Sigma_phi_phi=diag(((mu*(mu_star.k1-mu_star.k0)+(mu_dagger.k1-mu_dagger.k0))^2)+(mu^2)*trigamma(a.k1)+((1-mu)^2)*trigamma(b.k1)-trigamma(phi.k1))

  Sigma_beta_beta=t(X)%*%Tb%*%(Q.inv^2)%*%(Phi^2)%*%C.k1%*%Sigma_mu_mu%*%Tb%*%X
  Sigma_beta_gamma=t(X)%*%Tb%*%(Q.inv^2)%*%Phi%*%C.k1%*%Sigma_mu_phi%*%Tg%*%Z
  Sigma_gamma_gamma=t(Z)%*%Tg%*%(Q.inv^2)%*%C.k1%*%Sigma_phi_phi%*%Tg%*%Z

  Sigma=rbind(cbind(Sigma_beta_beta,Sigma_beta_gamma),cbind(t(Sigma_beta_gamma),Sigma_gamma_gamma))

  V=n*solve(Lambda)%*%Sigma%*%t(solve(Lambda))
  #V=n*MASS::ginv(Lambda)%*%Sigma%*%t(MASS::ginv(Lambda))

  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma
  result$Cov=V/n

  result$K=t(X)%*%Tb%*%Q.inv%*%(Phi^2)%*%C.k0%*%Lambda_mu_mu%*%Tb%*%X
  result$Std.Error=suppressWarnings(c(sqrt(diag(V/n))))

  return(result)
}

#Hat matrix
hatvalues.LSMLE=function(object)
{
  mu_hat=object$fitted.values$mu.predict
  phi_hat=object$fitted.values$phi.predict
  y=object$y
  X=object$model$mean
  linkobj=set.link2(link.mu=object$link,link.phi=object$link.phi)
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)

  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  W.PHI=diag(x=phi_hat*V_star*((d.link.mu)^(-2)))
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)
  return(diag(H))
}

