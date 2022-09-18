#' @export
SMLE.fit=function(y,x,z,alpha=NULL,link="logit",link.phi="log",control=robustbetareg.control(...),...)
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
    return(Opt.Tuning.SMLE(y,x,z,link,link.phi,control))
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
  theta=tryCatch(optim(par=start_theta,fn=L_q,gr=Psi_SMLE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1)),error=function(e){
    theta$msg<-e$message;check<<-F;return(theta)})
  # gc()
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
  M.SMLE=tryCatch(SMLE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj),error=function(e) NULL)
  if(is.null(M.SMLE))
  {
    vcov=std.error.SMLE=NaN
  }else{
    vcov=M.SMLE$Cov
    std.error.SMLE=M.SMLE$Std.Error
  }

  #Register of output values
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients#Coefficients Regression
  result$vcov=vcov#Expected Covariance Matrix
  pseudor2 = if (var(eta) * var(qlogis(y)) <= 0){NA}
  else{cor(eta, linkobj$linkfun.mu$linkfun(y))^2}
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.SMLE)))#Convergence
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)#Fitted Values
  result$start=start_theta #Started Point
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
  result$method="SMLE"
  if(any(is.na(std.error.SMLE)))
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
  if(!any(is.na(std.error.SMLE)))
  {#Standard Error
    names(std.error.SMLE)<-c(colnames(x),colnames(z))
    se.beta=std.error.SMLE[1:k]
    se.gamma=std.error.SMLE[1:m+k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error=coef.std.error
  }
  class(result)="robustbetareg"
  return(result)
}



###########################################################################################################################

# Auto Selecting tuning parameter algorithm
Opt.Tuning.SMLE=function(y,x,z,link,link.phi,control)
{
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  SMLE.list=SMLE.par=list()
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
    # est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi), control=betareg.control(start=Initial.points(y,x,z))),error=function(e) NULL)
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
    control$start=Est.param
    SMLE.par=tryCatch(SMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {SMLE.par$converged<-FALSE; return(SMLE.par)})
    if(!SMLE.par$converged || is.null(SMLE.par) || any(is.na(do.call("c",SMLE.par$coefficients)/do.call("c",SMLE.par$std.error))) || is.null(do.call("c",SMLE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    SMLE.list[[k]]<-SMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",SMLE.par$coefficients)/do.call("c",SMLE.par$std.error)))
  }
  sqv=as.numeric(SQV(zq.t,n,p))
  if(unstable || sqv==0)#Lack of stability
  {
    SMLE.par.star=SMLE.list[[1]]
    SMLE.par.star$sqv=sqv
    SMLE.par.star$Optimal.Tuning=TRUE
    SMLE.par.star$message="Lack of stability"
    return(SMLE.par.star)
  }
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    SMLE.par.star<-SMLE.list[[1]]
    SMLE.par.star$sqv=sqv
    SMLE.par.star$Optimal.Tuning=TRUE
    rm(SMLE.list)
    return(SMLE.par.star)
  }
  if(alpha.ind<8){#Which Tuning satisfy the condition os stability
    SMLE.par.star<-SMLE.list[[alpha.ind+1]]
    SMLE.par.star$sqv=sqv
    SMLE.par.star$Optimal.Tuning=TRUE
    rm(SMLE.list)
    return(SMLE.par.star)
  }

  reached=FALSE
  k=M1+1
  while(sqv.unstable & !reached)#Seek within the next grid of tuning
  {
    control$start=Est.param
    SMLE.par=tryCatch(SMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){SMLE.par$converged<-FALSE; return(SMLE.par)})
    if(!SMLE.par$converged || any(is.na(do.call("c",SMLE.par$coefficients)/do.call("c",SMLE.par$std.error))) || is.null(do.call("c",SMLE.par$std.error)) || !SMLE.par$converged )
    {
      unstable=T
      break
    }
    SMLE.list[[k]]=SMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",SMLE.par$coefficients)/do.call("c",SMLE.par$std.error)))
    sqv=as.numeric(SQV(zq.t,n,p))
    sqv.test=sqv[(k-M):(k-1)]
    if(all(sqv.test<=L))
    {
      sqv.unstable=F
    }
    k=k+1
    if(k>=K)#Condition of Step-6
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
    SMLE.par.star=SMLE.list[[1]]
    SMLE.par.star$sqv=sqv
    SMLE.par.star$Optimal.Tuning=TRUE
    SMLE.par.star$message="Lack of stability"
  }else{
    SMLE.par.star=SMLE.list[[(k-1-M)]]
    SMLE.par.star$sqv=sqv
    SMLE.par.star$Optimal.Tuning=TRUE
  }
  return(SMLE.par.star)
}


#Objective Function - SMLE Log-likelihood
L_q=function(theta,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]

  mu_q = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_q = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  phi <-(1/q)*(phi_q - 2) + 2
  mu <- ((1/q)*(mu_q*phi_q - 1) + 1)/phi
  a <- mu*phi
  b <- (1 - mu)*phi
  log_likS <- sum(dbeta(y, a, b, log = F)^(alpha))#function to be maximized
  return(log_likS)
}

#Gradiente: Psi_beta SMLE
Psi_Beta_SMLE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)

  phi_q=q*(phi_hat-2)+2
  mu_q=(q*(mu_hat*phi_hat-1)+1)/phi_q
  y_star=log(y)-log(1-y)
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  Tb = (link.model$link.mu$d.linkfun(mu_q))^(-1)

  result=t(X)%*%diag(phi_q*Tb/q)%*%(y_star-mu_star)
  return(result)
}

#Gradiente: Psi_gamma SMLE
Psi_Gamma_SMLE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)

  phi_q=q*(phi_hat-2)+2
  mu_q=(q*(mu_hat*phi_hat-1)+1)/phi_q
  y_dagger=log(1-y)
  y_star=log(y)-y_dagger

  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  mu_dagger=digamma((1-mu_hat)*phi_hat)-digmma(phi_hat)

  Tg = (link.model$linkfun.phi$d.linkfun(phi_q))^1

  result=t(Z)%*%diag(mu_q*Tg/q)%*%(mu_q*(y_star-mu_star)+(y_dagger-mu_dagger))
  return(result)
}

#Gradiente: Psi SMLE
Psi_SMLE=function(theta,y,X,Z,alpha,link_mu,link_phi){
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]
  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  psi_beta = Psi_Beta_SMLE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)
  psi_gamma = Psi_Gamma_SMLE(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi)

  return(c(psi_beta,psi_gamma))
}

#Sandwich Matrix - SMLE
SMLE_Cov_Matrix = function(muhat_q,phihat_q,X,Z,alpha,linkobj) {
  q=1-alpha
  ahat_q <- muhat_q*phihat_q
  bhat_q <- (1.0 - muhat_q)*phihat_q
  #mean link function
  T_1_q=diag((linkobj$linkfun.mu$d.linkfun(muhat_q))^(-1))
  #precision link funtion
  T_2_q=diag((linkobj$linkfun.phi$d.linkfun(phihat_q))^(-1))
  #Reparameterization
  phihat_n <-   (1.0/q)*(phihat_q - 2.0) + 2.0
  muhat_n <- 	((1.0/q)*(muhat_q*phihat_q  - 1.0) + 1.0)/phihat_n

  ahat_n <- muhat_n*phihat_n
  bhat_n <- (1.0 - muhat_n)*phihat_n
  phihat_2_q <- (2.0 - q)*(phihat_n - 2.0) + 2.0;	#expression of phi_(2-q)
  muhat_2_q <- ((2.0 - q)*(ahat_n - 1.0) + 1.0)/phihat_2_q; #expression of mu_(2-q)
  a2_qhat <- muhat_2_q*phihat_2_q
  b2_qhat <- (1.0 - muhat_2_q)*phihat_2_q

  mustarhat_n <- digamma(ahat_n) - digamma(bhat_n)
  mudaggerhat_n <-  digamma(bhat_n) - digamma(phihat_n)
  mustarhat_2_q <- digamma(a2_qhat) - digamma(b2_qhat)
  mudaggerhat_2_q <-  digamma(b2_qhat) - digamma(phihat_2_q)

  muhat_d_2_q <- muhat_q*(mustarhat_2_q - mustarhat_n) + mudaggerhat_2_q - mudaggerhat_n
  m_phiq <- diag(phihat_q)

  psi1_n <- trigamma(ahat_n)
  psi2_n <- trigamma(bhat_n)
  psi3_n <- trigamma(phihat_n)

  V_n <- diag(psi1_n + psi2_n)
  B1 <- diag(exp(q*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) - (lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))
  B2 <- diag(exp(lgamma(a2_qhat) + lgamma(b2_qhat) - lgamma(phihat_2_q) - (2.0*(1.0-q)*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n))  +
                                                                             lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))
  C_q_0 <- diag(phihat_q*(muhat_q*psi1_n - (1.0 - muhat_q)*psi2_n))
  D_q_0 <- diag((muhat_q^2.0)*psi1_n + ((1.0 - muhat_q)^2.0)*psi2_n - psi3_n)

  psi1_2_q <- trigamma(a2_qhat)
  psi2_2_q <- trigamma(b2_qhat)
  psi3_2_q <- trigamma(phihat_2_q)

  V_2_q <- diag(psi1_2_q + psi2_2_q)
  C_q_2_q <- diag(phihat_q*(muhat_q*psi1_2_q - (1.0 - muhat_q)*psi2_2_q))
  D_q_2_q <- diag((muhat_q^2.0)*psi1_2_q + ((1.0 - muhat_q)^2.0)*psi2_2_q - psi3_2_q)
  M1 <- diag(mustarhat_2_q - mustarhat_n)
  M2 <- diag(muhat_d_2_q)
  M3 <- diag((mustarhat_2_q - mustarhat_n)*muhat_d_2_q)

  #Matrix Jq
  Jq_betabeta <- as.matrix(t(X)%*%B1%*%(T_1_q^2.0)%*%(m_phiq^2.0)%*%V_n%*%X)
  Jq_betagamma <- as.matrix(t(X)%*%B1%*%T_1_q%*%T_2_q%*%C_q_0%*%Z)
  Jq_gammagamma <- as.matrix(t(Z)%*%B1%*%(T_2_q^2.0)%*%D_q_0%*%Z)
  Jq=-(q^(-1.0))*rbind(cbind(Jq_betabeta,Jq_betagamma),cbind(t(Jq_betagamma),Jq_gammagamma))

  #Matrix Kq
  Kq_betabeta <- as.matrix(t(X)%*%B2%*%(T_1_q^2.0)%*%(m_phiq^2.0)%*%(V_2_q+M1^2.0)%*%X)
  Kq_betagamma <- as.matrix(t(X)%*%B2%*%T_1_q%*%T_2_q%*%m_phiq%*%((C_q_2_q*(1/phihat_q)) + M3)%*%Z)
  Kq_gammagamma <- as.matrix(t(Z)%*%B2%*%(T_2_q^2.0)%*%(D_q_2_q+M2^2.0)%*%Z)
  Kq=(q^(-2.0))*rbind(cbind(Kq_betabeta,Kq_betagamma),cbind(t(Kq_betagamma),Kq_gammagamma))

  #Covariance Matrix
  Vq <- tryCatch( solve(Jq)%*%Kq%*%t(solve(Jq)), error=function(e) {e})  #asymptotic covariance matrix
  if(is.error(Vq)){
    Vq <- Ginv(Jq)%*%Kq%*%t(Ginv(Jq))
  }

  result=list()
  result$Lambda=Jq
  result$Sigma=Kq
  result$Cov=Vq
  result$Std.Error=suppressWarnings(t(sqrt(diag(Vq))))

  return(result)
}#ends Covariance matrix function

