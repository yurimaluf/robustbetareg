#' @rdname robustbetareg
#'
#' @usage MDPDE.fit(y,x,z, alpha=NULL,link="logit",
#' link.phi="log",control=robustbetareg.control(...),...)
#'
#' @param y,x,z \code{y} should be a numeric response vector (\eqn{y\in(0,1)}), \code{x} should be a numeric regresor matrix and \code{z} should be
#' a regressor matrix for the precision model, where there is the intercept only.
#' @param ... argument to be passed to \code{\link[=robustbetareg.control]{robustbetareg.control}}
#'
#' @importFrom stats dbeta
#' @importFrom stats cor
#' @importFrom stats var
#'
#'@export
MDPDE.fit=function(y,x,z,alpha=NULL,link="logit",link.phi="log",control=robustbetareg.control(...),...)
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
    return(Opt.Tuning.MDPDE(y,x,z,link,link.phi,control))
  }
  if(is.null(start_theta))
  {
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
    theta$x=start_theta
    theta$converged=mle$converged
  }
  #Point Estimation
  q=1-alpha
  check=TRUE
  theta=tryCatch(optim(par=start_theta,fn=D_q,gr=Psi_MDPDE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1)),error=function(e){
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
  M.MDPDE=tryCatch(MDPDE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj),error=function(e) NULL)
  if(is.null(M.MDPDE))
  {
    vcov=std.error.MDPDE=NaN
  }else{
    vcov=M.MDPDE$Cov
    std.error.MDPDE=M.MDPDE$Std.Error
  }

  #Register of output values
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients#Coefficients Regression
  result$vcov=vcov#Expected Covariance Matrix
  pseudor2 = if (var(eta) * var(qlogis(y)) <= 0){NA}
  else{cor(eta, linkobj$linkfun.mu$linkfun(y))^2}
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.MDPDE)))#Convergence
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)#Fitted Values
  result$start=start_theta#Started Point
  result$weights=(dbeta(y,mu_hat*phi_hat,(1-mu_hat)*phi_hat))^(alpha)#Weights
  result$Tuning=alpha#Tuning
  result$residuals=sweighted2_res(mu_hat,phi_hat,y=y,X=x,linkobj = linkobj)#Standard Residual Report
  result$n=length(mu_hat)#Sample Size
  result$link=link#Mean Link Function
  result$link.phi=link.phi#Precision Link Function
  result$Optimal.Tuning=alpha.optimal#Flag - Data-driven algorithm tuning selector
  result$pseudo.r.squared=pseudor2
  result$control=control#Extra options package
  result$call=match.call()
  result$method="MDPDE"
  if(any(is.na(std.error.MDPDE)))
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
  if(!any(is.na(std.error.MDPDE)))
  {#Standard Error
    names(std.error.MDPDE)<-c(colnames(x),colnames(z))
    se.beta=std.error.MDPDE[1:k]
    se.gamma=std.error.MDPDE[1:m+k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error=coef.std.error
  }
  class(result)="robustbetareg"
  return(result)
}


###########################################################################################################################

# Auto Selecting tuning parameter algorithm
Opt.Tuning.MDPDE=function(y,x,z,link,link.phi,control)
{
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  MDPDE.list=MDPDE.par=list()
  zq.t=NULL
  alpha_tuning=seq(0,0.3,0.02)
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
  p=length(Est.param)
  control$start=Est.param
  for(k in 1:M1)#First M1 attempts of tuning parameters
  {
    MDPDE.par=tryCatch(MDPDE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {MDPDE.par$converged<-FALSE; return(MDPDE.par)})
    if(!MDPDE.par$converged || is.null(MDPDE.par) || any(is.na(do.call("c",MDPDE.par$coefficients)/do.call("c",MDPDE.par$std.error))) || is.null(do.call("c",MDPDE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    MDPDE.list[[k]]<-MDPDE.par
    zq.t=unname(rbind(zq.t,do.call("c",MDPDE.par$coefficients)/do.call("c",MDPDE.par$std.error)))
  }

  sqv=as.numeric(SQV(zq.t,n,p))
  if(unstable || sqv==0)#Lack of stability
  {
    MDPDE.par.star=MDPDE.list[[1]]
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
    MDPDE.par.star$message="Lack of stability"
    return(MDPDE.par.star)
  }
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    MDPDE.par.star<-MDPDE.list[[1]]
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
    rm(MDPDE.list)
    return(MDPDE.par.star)
  }

  if(alpha.ind<8){#Which Tuning satisfy the condition os stability
    MDPDE.par.star<-MDPDE.list[[alpha.ind+1]]
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
    rm(MDPDE.list)
    return(MDPDE.par.star)
  }

  reached=FALSE
  k=M1+1
  while(sqv.unstable & !reached)#Seek within the next grid of tuning
  {
    MDPDE.par=tryCatch(MDPDE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){MDPDE.par$converged<-FALSE; return(MDPDE.par)})
    if(!MDPDE.par$converged || any(is.na(do.call("c",MDPDE.par$coefficients)/do.call("c",MDPDE.par$std.error))) || is.null(do.call("c",MDPDE.par$std.error)) || !MDPDE.par$converged )
    {
      unstable=T
      break
    }
    MDPDE.list[[k]]=MDPDE.par
    zq.t=unname(rbind(zq.t,do.call("c",MDPDE.par$coefficients)/do.call("c",MDPDE.par$std.error)))
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
    k=suppressWarnings(max(1,min(which(zoo::rollapply(sqv<L,M,sum)==M)))+M+1)
  }
  if(k>=K || unstable)
  {
    MDPDE.par.star=MDPDE.list[[1]]
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
    MDPDE.par.star$message="Lack of stability"
  }else{
    MDPDE.par.star=MDPDE.list[[(k-1-M)]]
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
  }
  return(MDPDE.par.star)
}


#Objective Function - MDPDE Log-likelihood
D_q=function(theta,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]

  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  if(alpha==0){
    D_q=sum(dbeta(y,mu_hat*phi_hat,(1-mu_hat)*phi_hat,log = T))
  }else{
    a0 = mu_hat*phi_hat
    b0 = (1-mu_hat)*phi_hat
    a_alpha = a0*(1+alpha)-alpha
    b_alpha = b0*(1+alpha)-alpha
    E_alpha =  exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) -  (1 + alpha)* (lgamma(a0) + lgamma(b0) - lgamma(a0+b0)))
    D_q = sum((1+alpha)/(alpha)*dbeta(y,mu_hat*phi_hat,(1-mu_hat)*phi_hat)^(alpha)-E_alpha)
  }
  return(D_q)
}

#Gradiente: Psi_beta MDPDE
Psi_Beta_MDPDE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)-alpha
  b_alpha=b0*(1+alpha)-alpha
  y_star = log(y)-log(1-y)
  Phi_Tb = diag(phi_hat*(link.model$linkfun.mu$d.linkfun(mu_hat))^(-1))

  mu_star = digamma(a0)-digamma(b0)
  mu_star_alpha = digamma(a_alpha)-digamma(b_alpha)
  Ubeta = y_star-mu_star
  E_Ubeta = mu_star_alpha-mu_star

  Walpha = diag(dbeta(y,mu_hat*phi_hat,(1-mu_hat)*phi_hat)^(alpha))
  Calpha =  exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) -  (1 + alpha)* (lgamma(a0) + lgamma(b0) - lgamma(a0+b0)))

  result = (1+alpha)*t(X)%*%Phi_Tb%*%(Walpha%*%Ubeta-Calpha%*%E_Ubeta)
  return(result)
}


#Gradiente: Psi_beta MDPDE
Psi_Gamma_MDPDE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)-alpha
  b_alpha=b0*(1+alpha)-alpha
  y_dagger = log(1-y)
  y_star = log(y)-y_dagger
  Tg = diag((link.model$linkfun.phi$d.linkfun(phi_hat))^(-1))

  mu_star = digamma(a0)-digamma(b0)
  mu_star_alpha = digamma(a_alpha)-digamma(b_alpha)
  mu_dagger = digamma(b0)-digamma(phi_hat)
  mu_dagger_alpha = digamma(b_alpha)-digamma(a_alpha+b_alpha)

  Ugamma = mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger)
  E_Ugamma = mu_hat*(mu_star_alpha-mu_star)+(mu_dagger_alpha-mu_dagger)

  Walpha = diag(dbeta(y,mu_hat*phi_hat,(1-mu_hat)*phi_hat)^(alpha))
  Calpha =  exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) -  (1 + alpha)* (lgamma(a0) + lgamma(b0) - lgamma(a0+b0)))

  result = (1+alpha)*t(Z)%*%Tg%*%(Walpha%*%Ugamma-Calpha%*%E_Ugamma)
  return(result)

}


#Gradiente: Psi MDPDE
Psi_MDPDE=function(theta,y,X,Z,alpha,link_mu,link_phi){
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]
  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  psi_beta = Psi_Beta_MDPDE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)
  psi_gamma = Psi_Gamma_MDPDE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)

  return(c(psi_beta,psi_gamma))
}


#Sandwich Matrix - MDPDE
MDPDE_Cov_Matrix = function(mu,phi,X,Z,alpha,linkobj){
  q_const=1-alpha
  a <- mu*phi
  b <- (1.0 - mu)*phi
  t_1 <- ((linkobj$linkfun.mu$d.linkfun(mu))^(-1))
  t_2 <- ((linkobj$linkfun.phi$d.linkfun(phi))^(-1))
  mustar <- psigamma(a, 0) - psigamma(b, 0)
  mudagger <- psigamma(b, 0) - psigamma(a + b, 0)
  m_phi <- diag(phi)
  psi1 <- psigamma(a, 1.0)
  psi2 <- psigamma(b, 1.0)
  psi3 <- psigamma(a + b, 1.0)
  a_q <- (2.0 - q_const)*(a - 1.0) + 1.0
  b_q <- (2.0 - q_const)*(b - 1.0) + 1.0
  psi1_q <- psigamma(a_q, 1.0)
  psi2_q <- psigamma(b_q, 1.0)
  psi3_q <- psigamma(a_q + b_q, 1.0)
  a2_q <- (3.0 - 2.0*q_const)*(a - 1.0) + 1.0
  b2_q <- (3.0 - 2.0*q_const)*(b - 1.0) + 1.0
  psi1_2q <- psigamma(a2_q, 1.0)
  psi2_2q <- psigamma(b2_q, 1.0)
  psi3_2q <- psigamma(a2_q + b2_q, 1.0)
  mustar_q <- psigamma(a_q, 0) - psigamma(b_q, 0)
  mustar_2q <- psigamma(a2_q, 0) - psigamma(b2_q, 0)
  mudagger_q <-  psigamma(b_q, 0) - psigamma(a_q + b_q, 0)
  mudagger_2q <-  psigamma(b2_q, 0) - psigamma(a2_q + b2_q, 0)
  K <-  exp(lgamma(a_q) + lgamma(b_q) - lgamma(a_q + b_q) - (2.0 - q_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))
  K2 <-  exp(lgamma(a2_q) + lgamma(b2_q) - lgamma(a2_q + b2_q) - (3.0 - 2.0*q_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))
  gama11_q <- diag(K*(phi^2.0)*(t_1^2.0)*(psi1_q + psi2_q + (mustar_q - mustar)^2.0))
  gama12_q <- diag(K*phi*t_1*t_2*(mu*(psi1_q +  psi2_q + (mustar_q - mustar)^2.0) - psi2_q + (mustar_q - mustar)*(mudagger_q - mudagger)))
  gama22_q <- diag(K*(t_2^2.0)*((mu^2.0)*psi1_q + ((1.0 - mu)^2.0)*psi2_q - psi3_q + (mu*(mustar_q - mustar) + mudagger_q - mudagger)^2.0))
  gama11_2q <- diag(K2*(phi^2.0)*(t_1^2.0)*(psi1_2q + psi2_2q + (mustar_2q - mustar)^2.0))
  gama12_2q <- diag(K2*phi*t_1*t_2*(mu*(psi1_2q + psi2_2q + (mustar_2q - mustar)^2.0) - psi2_2q + (mustar_2q - mustar)*(mudagger_2q - mudagger)))
  gama22_2q <- diag(K2*(t_2^2.0)*((mu^2.0)*psi1_2q + ((1.0 - mu)^2.0)*psi2_2q - psi3_2q + (mu*(mustar_2q - mustar) + mudagger_2q - mudagger)^2.0))
  E1q <- diag(K*phi*t_1*(mustar_q - mustar))
  E2q <- diag(K*t_2*(mu*(mustar_q - mustar) + mudagger_q - mudagger))

  #Matrix Psin
  Psin_betabeta <- -(2.0 - q_const)*as.matrix(t(X)%*%gama11_q%*%X)
  Psin_betagamma <- -(2.0 - q_const)*as.matrix(t(X)%*%gama12_q%*%Z)
  Psin_gammagamma <- -(2.0 - q_const)*as.matrix(t(Z)%*%gama22_q%*%Z)
  Psin=rbind(cbind(Psin_betabeta,Psin_betagamma),cbind(t(Psin_betagamma),Psin_gammagamma))

  #Matrix Omegan
  Omegan_betabeta <- ((2.0 - q_const)^2.0)*as.matrix(t(X)%*%(gama11_2q - E1q^2.0)%*%X)
  Omegan_betagamma <- ((2.0 - q_const)^2.0)*as.matrix(t(X)%*%(gama12_2q - E1q*E2q)%*%Z)
  Omegan_gammagamma <- ((2.0 - q_const)^2.0)*as.matrix(t(Z)%*%(gama22_2q - E2q^2.0)%*%Z)
  Omegan=rbind(cbind(Omegan_betabeta,Omegan_betagamma),cbind(t(Omegan_betagamma),Omegan_gammagamma))

  inv.Psin=tryCatch(solve(Psin), error=function(e) {e})
  if(!BBmisc::is.error(inv.Psin)){
    Vq <- inv.Psin%*%(Omegan)%*%t(inv.Psin)
  }else{
    inv.Psin=MASS::ginv(Psin)
    Vq <- inv.Psin%*%(Omegan)%*%t(inv.Psin)
  }

  result=list()
  result$Lambda=Psin
  result$Sigma=Omegan
  result$Cov=Vq
  result$Std.Error=suppressWarnings(t(sqrt(diag(Vq))))

  return(result)
}#ends Covariance matrix function

