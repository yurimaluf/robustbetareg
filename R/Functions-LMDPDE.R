#' @export
LMDPDE.fit=function(y,x,z,alpha=NULL,link="logit",link.phi="log",control=robustbetareg.control(...), ...)
{
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
    return(Opt.Tuning.LMDPDE(y,x,z,link,link.phi,control))
  }
  if(is.null(start_theta))
  {
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
  }
  #Point Estimation
  check=TRUE
  theta=tryCatch(optim(par=start_theta,fn=D_alpha_R,gr=Psi_LMDPDE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1,maxit=10000)),error=function(e){
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
  MM=tryCatch(LMDPDE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj),error=function(e) NULL)
  if(is.null(MM))
  {
    vcov= std.error.LMDPDE=NULL
  }else{
    vcov=MM$Cov
    std.error.LMDPDE=MM$Std.Error
  }
  #Register of output values
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients#Coefficients Regression
  result$vcov=vcov#Expected Covariance Matrix
  pseudor2 <- if (var(eta) * var(qlogis(y)) <= 0){NA}
  else{cor(eta, linkobj$linkfun.mu$linkfun(y))^2}
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.LMDPDE)))#Convergence
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)#Fitted Values
  result$start=start_theta#Started Point
  result$weights=degbeta(y_star,mu_hat,phi_hat)^(alpha)#Weights
  result$Tuning=alpha#Tuning
  result$residuals=sweighted2_res(mu_hat,phi_hat,y=y,X=x,linkobj = linkobj)#Standard Residual Report
  result$n=length(mu_hat)#Sample Size
  result$link=link#Mean Link Function
  result$link.phi=link.phi#Precision Link Function
  result$Optimal.Tuning=alpha.optimal#Flag - Data-driven algorithm tuning selector
  result$pseudo.r.squared=pseudor2
  result$control=control#Extra options package
  result$call=match.call()
  result$method="LMDPDE"
  if(any(is.na(std.error.LMDPDE)))
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
  if(!any(is.na(std.error.LMDPDE)))
  {#Standard Error
    names(std.error.LMDPDE)<-c(colnames(x),colnames(z))
    se.beta=std.error.LMDPDE[1:k]
    se.gamma=std.error.LMDPDE[1:m+k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error=coef.std.error
  }
  class(result)="robustbetareg"
  return(result)
}


###########################################################################################################

# Auto Selecting tuning parameter algorithm
Opt.Tuning.LMDPDE=function(y,x,z,link,link.phi,control)
{
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LMDPDE.list=LMDPDE.par=list()
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
    # est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi), control=betareg.control(start=Initial.points(y,x,z))),error=function(e) NULL)
  }
  if(!is.null(est.log.lik)){
    Est.param=do.call("c",est.log.lik$coefficients)
    names(Est.param)=c(colnames(x),colnames(z))
  }else{
    Est.param=Initial.points(y,x,z)
  }
  ponto.inicial.robst=Initial.points(y,x,z)
  ponto.inicial.temp=Est.param
  names(ponto.inicial.robst)=names(ponto.inicial.temp)=c(colnames(x),colnames(z))
  p=length(Est.param)
  for(k in 1:M1)#First M1 attempts of tuning parameters
  {
    control$start=ponto.inicial.temp
    LMDPDE.par=tryCatch(LMDPDE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(LMDPDE.par$converged)
    {
      ponto.inicial.temp=c(LMDPDE.par$coefficients$mean,LMDPDE.par$coefficients$precision)
    }
    if(!LMDPDE.par$converged || is.null(LMDPDE.par) || any(is.na(do.call("c",LMDPDE.par$coefficients)/do.call("c",LMDPDE.par$std.error))) || is.null(do.call("c",LMDPDE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LMDPDE.list[[k]]<-LMDPDE.par
    zq.t=unname(rbind(zq.t,do.call("c",LMDPDE.par$coefficients)/do.call("c",LMDPDE.par$std.error)))
  }
  #sqv=as.numeric(SQV_Cpp(zq.t,n,p))
  sqv=as.numeric(SQV(zq.t,n,p))
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    LMDPDE.par.star<-LMDPDE.list[[1]]
    LMDPDE.par.star$sqv=sqv
    LMDPDE.par.star$Optimal.Tuning=TRUE
    rm(LMDPDE.list)
    return(LMDPDE.par.star)
  }
  if(unstable)#Lack of stability
  {
    LMDPDE.par.star=LMDPDE.list[[1]]
    LMDPDE.par.star$sqv=sqv
    LMDPDE.par.star$Optimal.Tuning=TRUE
    LMDPDE.par.star$message="Lack of stability"
    return(LMDPDE.par.star)
  }
  if(alpha.ind<8){#Which Tuning satisfy the condition os stability
    LMDPDE.par.star<-LMDPDE.list[[alpha.ind+1]]
    LMDPDE.par.star$sqv=sqv
    LMDPDE.par.star$Optimal.Tuning=TRUE
    rm(LMDPDE.list)
    return(LMDPDE.par.star)
  }

  reached=FALSE
  k=M1+1
  while(sqv.unstable & !reached)#Seek within the next grid of tuning
  {
    LMDPDE.par=tryCatch(LMDPDE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(LMDPDE.par$converged)
    {
      ponto.inicial.temp=c(LMDPDE.par$coefficients$mean,LMDPDE.par$coefficients$precision)
    }
    if(any(is.na(do.call("c",LMDPDE.par$coefficients)/do.call("c",LMDPDE.par$std.error))) || is.null(do.call("c",LMDPDE.par$std.error)))
    {
      unstable=T
      break
    }
    LMDPDE.list[[k]]=LMDPDE.par
    zq.t=unname(rbind(zq.t,do.call("c",LMDPDE.par$coefficients)/do.call("c",LMDPDE.par$std.error)))
    #sqv=as.numeric(SQV_Cpp(zq.t,n,p))
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
    LMDPDE.par.star=LMDPDE.list[[1]]
    LMDPDE.par.star$sqv=sqv
    LMDPDE.par.star$Optimal.Tuning=TRUE
    LMDPDE.par.star$message="Lack of stability"
  }else{
    LMDPDE.par.star=LMDPDE.list[[(k-1-M)]]
    LMDPDE.par.star$sqv=sqv
    LMDPDE.par.star$Optimal.Tuning=TRUE
  }
  return(LMDPDE.par.star)
}


#Objective Function - LMDPDE Log-likelihood
D_alpha_R=function(theta,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  k=ncol(X)
  m=ncol(Z)
  y_star=log(y)-log(1-y)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]
  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  if(alpha==0){
    D_q = sum(log(degbeta(y_star,mu_hat,phi_hat)))
  }else{
    a0 = mu_hat*phi_hat
    b0 = (1-mu_hat)*phi_hat
    a_alpha = a0*(1+alpha)
    b_alpha = b0*(1+alpha)
    E_alpha = exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0))
    D_q = sum((1+alpha)/(alpha)*degbeta(y_star,mu_hat,phi_hat)^(alpha)-E_alpha)
  }
  return(D_q)
}

#Gradiente: Psi_beta LMDPDE
Psi_Beta_LMDPDE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)
  b_alpha=b0*(1+alpha)
  y_star=log(y)-log(1-y)
  Phi_Tb = phi_hat*(link.model$linkfun.mu$d.linkfun(mu_hat))^(-1)

  mu_star = digamma(a0)-digamma(b0)
  mu_star_alpha = digamma(a_alpha)-digamma(b_alpha)
  Ubeta = y_star-mu_star
  E_Ubeta = mu_star_alpha-mu_star

  Walpha = diag(Phi_Tb*degbeta(y_star,mu_hat,phi_hat)^(alpha))
  Calpha = diag(Phi_Tb*exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0)))

  result = (1+alpha)*t(X)%*%Phi_Tb%*%(Walpha%*%Ubeta-Calpha%*%E_Ubeta)
  return(result)
}

#Gradiente: Psi_gamma LMDPDE
Psi_Gamma_LMDPDE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)
  b_alpha=b0*(1+alpha)
  y_dagger = log(1-y)
  y_star = log(y)-y_dagger
  Tg = diag((link.model$linkfun.phi$d.linkfun(phi_hat))^(-1))

  mu_star = digamma(a0)-digamma(b0)
  mu_star_alpha = digamma(a_alpha)-digamma(b_alpha)
  mu_dagger = digamma(b0)-digamma(phi_hat)
  mu_dagger_alpha = digamma(b_alpha)-digamma(phi_hat*(1+alpha))
  Ugamma = mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger)
  E_Ugamma = mu_hat*(mu_star_alpha-mu_star)+(mu_dagger_alpha-mu_dagger)

  Walpha = degbeta(y_star,mu_hat,phi_hat)^(alpha)
  Calpha = exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lebeta(a0,b0))

  result = (1+alpha)*t(Z)%*%Tg%*%(Walpha%*%Ugamma-Calpha%*%E_Ugamma)
  return(result)
}

#Gradiente: Psi LMDPDE
Psi_LMDPDE=function(theta,y,X,Z,alpha,link_mu,link_phi){
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]
  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  psi_beta = Psi_Beta_LMDPDE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)
  psi_gamma = Psi_Gamma_LMDPDE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)

  return(c(psi_beta,psi_gamma))
}


# Sandwich Matrix - LMDPDE
LMDPDE_Cov_Matrix=function(mu,phi,X,Z,alpha,linkobj)
{
  n=length(mu)
  a0=mu*phi
  b0=(1-mu)*phi
  a_alpha=(1+alpha)*a0
  b_alpha=(1+alpha)*b0
  a_2alpha=(1+2*alpha)*a0
  b_2alpha=(1+2*alpha)*b0

  K_alpha=exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0))
  K_2alpha=exp(lbeta(a_2alpha,b_2alpha)-(1+2*alpha)*lbeta(a0,b0))

  #mean link function
  Tb=diag((linkobj$linkfun.mu$d.linkfun(mu))^(-1))
  #precision link funtion
  Tg=diag((linkobj$linkfun.phi$d.linkfun(phi))^(-1))

  mu_alpha_star=digamma(a_alpha)-digamma(b_alpha)
  mu_alpha_dagger=digamma(b_alpha)-digamma(a_alpha+b_alpha)
  mu_star=digamma(a0)-digamma(b0)
  mu_dagger=digamma(b0)-digamma(a0+b0)
  nu_alpha=mu^2*trigamma(a_alpha)+(1-mu)^2*trigamma(b_alpha)-trigamma(a_alpha+b_alpha)

  mu_2alpha_star=digamma(a_2alpha)-digamma(b_2alpha)
  mu_2alpha_dagger=digamma(b_2alpha)-digamma(a_2alpha+b_2alpha)
  kappa_mu=phi*K_alpha*(mu_alpha_star-mu_star)
  kappa_phi=K_alpha*(mu*(mu_alpha_star-mu_star)+mu_alpha_dagger-mu_dagger)
  E.u2.phi_2alpha=(mu*(mu_2alpha_star-mu_star)+(mu_2alpha_dagger-mu_dagger))^2+mu^2*trigamma(a_2alpha)+(1-mu)^2*trigamma(b_2alpha)-trigamma(a_2alpha+b_2alpha)

  Lambda_mu_mu=diag(phi^2*K_alpha*((mu_alpha_star-mu_star)^2+trigamma(a_alpha)+trigamma(b_alpha)))
  Lambda_mu_phi=diag(phi*K_alpha*(mu*trigamma(a_alpha)-(1-mu)*trigamma(b_alpha)+mu*(mu_alpha_star-mu_star)^2+(mu_alpha_star-mu_star)*(mu_alpha_dagger-mu_dagger)))
  Lambda_phi_phi=diag(K_alpha*((mu*(mu_alpha_star-mu_star)+mu_alpha_dagger-mu_dagger)^2+nu_alpha))

  Sigma_mu_mu=diag(phi^2*K_2alpha*(trigamma(a_2alpha)+trigamma(b_2alpha)+(mu_2alpha_star-mu_star)^2)-kappa_mu^2)
  Sigma_mu_phi=diag(K_2alpha*phi*(mu*trigamma(a_2alpha)-(1-mu)*trigamma(b_2alpha)+mu*(mu_2alpha_star-mu_star)^2+(mu_2alpha_star-mu_star)*(mu_2alpha_dagger-mu_dagger))-kappa_mu*kappa_phi)
  Sigma_phi_phi=diag(K_2alpha*E.u2.phi_2alpha-kappa_phi^2)

  Lambda_beta_beta=t(X)%*%Tb%*%Lambda_mu_mu%*%Tb%*%X
  Lambda_beta_gamma=t(X)%*%Tb%*%Lambda_mu_phi%*%Tg%*%Z
  Lambda_gamma_gamma=t(Z)%*%Tg%*%Lambda_phi_phi%*%Tg%*%Z

  Lambda=rbind(cbind(Lambda_beta_beta,Lambda_beta_gamma),cbind(t(Lambda_beta_gamma),Lambda_gamma_gamma))

  Sigma_beta_beta=t(X)%*%Tb%*%Sigma_mu_mu%*%Tb%*%X
  Sigma_beta_gamma=t(X)%*%Tb%*%Sigma_mu_phi%*%Tg%*%Z
  Sigma_gamma_gamma=t(Z)%*%Tg%*%Sigma_phi_phi%*%Tg%*%Z

  Sigma=rbind(cbind(Sigma_beta_beta,Sigma_beta_gamma),cbind(t(Sigma_beta_gamma),Sigma_gamma_gamma))

  V=n*solve(Lambda)%*%Sigma%*%t(solve(Lambda))
  #V=n*MASS::ginv(Lambda)%*%Sigma%*%t(MASS::ginv(Lambda))

  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma

  result$Cov=V/n
  result$Std.Error=suppressWarnings(c(sqrt(diag(V/n))))

  return(result)
}

