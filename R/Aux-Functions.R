#' Create a Link for robust beta regression model
#'
#' This function provides several link functions for robust beta regression models
#'
#' @param link.mu character; one of "logit"(default), "probit", "cauchit", "cloglog" and "loglog"
#' @param link.phi character; one of "log"(default) "indentity" and "sqrt"
#'
#' @return A structure with link function, inverse link function, derivative \eqn{d\eta_{\mu}/d\mu} and \eqn{d^2\eta_{\mu}/d\mu^2} (\eqn{d\eta_{\phi}/d\phi} and \eqn{d^2\eta_{\phi}/d\phi^2}), for both models: mean and precision,
#' where \eqn{\eta_{\mu}=\textbf{X}^{\top}\beta} and \eqn{\eta_{\phi}=\textbf{Z}^{\top}\gamma}.
#'
#' @examples
#' \dontrun{
#' links=set.link(link.mu="cauchit",link.phi="sqrt")
#' attributes(links)}
#'
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom miscTools ddnorm
#'
set.link=function(link.mu="logit",link.phi="log")
{
  #Mean Links
  switch(link.mu,
         logit={
           linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(log(mu)-log(1-mu))
           }
           d.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return((mu-mu^2)^(-1))
           }
           d2.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return((2*mu-1)/(mu*(1-mu))^2)
           }
           inv.link=function(eta)
           {
             return(as.numeric(pmax(pmin(exp(eta-Rmpfr::log1pexp(eta)),1-.Machine$double.eps),.Machine$double.eps)))
           }
         },
         probit={
           linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(qnorm(mu))
           }
           d.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return((dnorm(qnorm(mu)))^(-1))
           }
           d2.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(-ddnorm(qnorm(mu))/(dnorm(qnorm(mu)))^3)
           }
           inv.link=function(eta)
           {
             return(as.numeric(pmax(pmin(pnorm(eta),1-.Machine$double.eps),.Machine$double.eps)))
           }
         },
         cloglog={
           linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(log(-log(1-mu)))
           }
           d.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return((-(1-mu)*log(1-mu))^(-1))
           }
           d2.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(-(log(1-mu)+1)/((1-mu)*log(1-mu))^2)
           }
           inv.link=function(eta)
           {
             return(as.numeric(pmax(pmin(1-exp(-exp(-eta)),1-.Machine$double.eps),.Machine$double.eps)))
           }
         },
         cauchit = {
           linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(tan(pi*(mu-0.5)))
           }
           d.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(pi*pracma::sec(pi*(mu-0.5))^2)
           }
           d2.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(2*pi*tan(pi*(mu-0.5))*pracma::sec(pi*(mu-0.5))^2)
           }
           inv.link=function(eta)
           {
             return(as.numeric(pmax(pmin(0.5+atan(eta)/pi,1-.Machine$double.eps),.Machine$double.eps)))
           }
         },
         loglog = {
           linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(-log(-log(mu)))
           }
           d.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return(-(mu*log(mu))^(-1))
           }
           d2.linkfun=function(mu)
           {
             mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
             return((log(mu)+1)/(mu*log(mu))^2)
           }
           inv.link=function(eta)
           {
             return(as.numeric(pmax(pmin(exp(-exp(-eta)),1-.Machine$double.eps),.Machine$double.eps)))
           }
         }, stop(gettextf("%s link.mu not recognised", sQuote(link.mu)), domain = NA))

  Linkfun.Mu=list(linkfun=linkfun,d.linkfun=d.linkfun,d2.linkfun=d2.linkfun,inv.link=inv.link)
  ### Precision Links
  switch(link.phi,
         log = {
           linkfun.phi=function(phi)
           {
             phi=pmax(phi,.Machine$double.eps)
             return(log(phi))
           }
           d.linkfun.phi=function(phi)
           {
             phi=pmax(phi,.Machine$double.eps)
             return((phi)^(-1))
           }
           d2.linkfun.phi=function(phi)
           {
             phi=pmax(phi,.Machine$double.eps)
             return(-(phi^(-2)))
           }
           inv.link.phi=function(eta)
           {
             return(pmax(as.numeric(exp(eta)),.Machine$double.eps))
           }
         },
         identity = {
           linkfun.phi=function(phi)
           {
             phi=pmax(phi,.Machine$double.eps)
             return(phi)
           }
           d.linkfun.phi=function(phi)
           {
             return(rep(1,length(phi)))
           }
           d2.linkfun.phi=function(phi)
           {
             return(rep(0,length(phi)))
           }
           inv.link.phi=function(eta)
           {
             return(as.numeric(eta))
           }
         },
         sqrt = {
           linkfun.phi=function(phi)
           {
             phi=pmax(phi,.Machine$double.eps)
             return(sqrt(phi))
           }
           d.linkfun.phi=function(phi)
           {
             return((2*sqrt(phi))^(-1))
           }
           d2.linkfun.phi=function(phi)
           {
             return(-0.25*(phi^(-3/2)))
           }
           inv.link.phi=function(eta)
           {
             return(pmax(as.numeric(eta^2),.Machine$double.eps))
           }
         }, stop(gettextf("%s link.phi not recognised", sQuote(link.phi)), domain = NA))

  Linkfun.Phi=list(linkfun=linkfun.phi,d.linkfun=d.linkfun.phi,d2.linkfun=d2.linkfun.phi,inv.link=inv.link.phi)
  ###
  linkobj=structure(list(linkfun.mu=Linkfun.Mu,linkfun.phi=Linkfun.Phi),name.link.mu=link.mu,name.link.phi=link.phi,class="link-rbr")
  return(linkobj)
}

#################################################################################################################################

#'
pearson_res=function(mu_hat,phi_hat,y)
{
  var.y=mu_hat*(1-mu_hat)/(1+phi_hat)
  ri=(y-mu_hat)/sqrt(var.y)
  return(ri)
}

#'
sweighted_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  nu=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  ri=(y_star-mu_star)/sqrt(nu)#standardized weighted residuals

  return(ri)
}

#'
sweighted2_res=function(mu_hat,phi_hat,y,X,linkobj)
{
  n=length(mu_hat)
  m=length(phi_hat)
  if(m==1){phi_hat=rep(phi_hat,n)}
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  y_star=log(y)-log(1-y)
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)

  W.PHI=diag(x=phi_hat*V_star*((d.link.mu)^(-2)))

  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)

  nu=V_star*(1-diag(H))
  diff=(y_star-mu_star)
  ri=diff/sqrt(nu)#standardized weighted residuals

  return(ri)
}

#'
sweighted3_res=function(mu_hat,phi_hat,alpha,y,X,linkobj)
{
  n=length(mu_hat)
  m=length(phi_hat)
  q=1-alpha
  if(m==1){phi_hat=rep(phi_hat,n)}
  a_0=mu_hat*phi_hat
  b_0=(1-mu_hat)*phi_hat
  a_alpha=a_0*(1+alpha)
  b_alpha=b_0*(1+alpha)

  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  y_star=log(y)-log(1-y)
  mu_star=digamma(a_alpha)-digamma(b_alpha)
  v_alpha=trigamma(a_alpha)+trigamma(b_alpha)
  c_alpha=(beta(a_alpha,b_alpha)^q)/beta(a_0,b_0)
  w_alpha=(degb(y_star,mu_hat,phi_hat/q))^(alpha)
  K_alpha=diag(v_alpha*c_alpha*(phi_hat/q)/d.link.mu^2)
  H_alpha=sqrt(K_alpha)%*%X%*%solve(t(X)%*%K_alpha%*%X)%*%t(X)%*%sqrt(K_alpha)

  nu=c_alpha*v_alpha*(1-diag(H_alpha))/q
  ri=(y_star-mu_star)*w_alpha/sqrt(nu)#standardized weighted residuals

  return(ri)
}

#'
weighted_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)

  nu=V_star*phi_hat
  ri=(y_star-mu_star)/sqrt(nu)#standardized weighted residuals

  return(ri)
}

#'
sweighted.gamma_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  xi=trigamma((1-mu_hat)*phi_hat)-trigamma(phi_hat)
  a=mu_hat*(y_star-mu_star)+log(1-y)-(digamma((1-mu_hat)*phi_hat)-digamma(phi_hat))
  ri=a/sqrt(xi) #standardized weighted residuals gamma
  return(ri)
}

#'
sweighted2.gamma_res=function(mu_hat,phi_hat,y,Z,linkobj)
{
  y_star=log(y/(1-y))
  n=length(mu_hat)
  m=length(phi_hat)
  if(m==1){phi_hat=rep(phi_hat,n)}
  d.link.phi=linkobj$linkfun.phi$d.linkfun(phi_hat)
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  xi=(trigamma(mu_hat*phi_hat)*mu_hat^2+trigamma((1-mu_hat)*phi_hat)*(1-mu_hat)^2-trigamma(phi_hat))
  Z=as.matrix(Z)
  D=diag(xi*(d.link.phi)^(-2))##Para funcao ligacao de phi
  G=sqrt(D)%*%Z%*%solve(t(Z)%*%D%*%Z)%*%t(Z)%*%sqrt(D)

  a=mu_hat*(y_star-mu_star)+log(1-y)-(digamma((1-mu_hat)*phi_hat)-digamma(phi_hat))
  nu=xi*(1-diag(G))
  ri=a/sqrt(nu) #standardized weighted residuals

  return(ri)
}

#'
combined_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  a=mu_hat*(y_star-mu_star)+log(1-y)-(digamma((1-mu_hat)*phi_hat)-digamma(phi_hat))
  zeta=trigamma(mu_hat*phi_hat)*(1+mu_hat)^2+trigamma((1-mu_hat)*phi_hat)*mu_hat^2-trigamma(phi_hat)
  ri=(y_star-mu_star+a)/sqrt(zeta) #combined residuals
  return(ri)
}

#'
combined.projection_res=function(mu_hat,phi_hat,y,X,Z,linkobj)
{
  n=length(mu_hat)
  m=length(phi_hat)
  if(m==1){phi_hat=rep(phi_hat,n)}
  #derivada da inversa da funcao ligacao
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  d.link.phi=linkobj$linkfun.phi$d.linkfun(phi_hat)

  y_star=log(y/(1-y))
  res.beta=sweighted_res(mu_hat,phi_hat,y)
  res.gamma=sweighted.gamma_res(mu_hat,phi_hat,y)
  v=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  b=mu_hat*v-trigamma((1-mu_hat)*phi_hat)
  xi=trigamma(mu_hat*phi_hat)*mu_hat^2+trigamma((1-mu_hat)*phi_hat)*(1-mu_hat)^2-trigamma(phi_hat)
  w=phi_hat*v*((d.link.mu)^(-2))
  TT=diag((d.link.mu)^(-1))
  HH=diag((d.link.phi)^(-1))
  B=diag(b)

  W.PHI=diag(x=phi_hat*w)
  W_PHI=diag(x=phi_hat/w)
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)

  D=diag(xi*(d.link.phi)^(-2))
  G=sqrt(D)%*%Z%*%solve(t(Z)%*%D%*%Z)%*%t(Z)%*%sqrt(D)
  M=sqrt(W_PHI)%*%TT%*%B%*%HH%*%sqrt(D)

  h=1-diag(H)
  g=1-diag(G)
  m=diag(M)
  nu=pmax(h+g+2*m,.Machine$double.eps)
  ri=(res.beta+res.gamma)/sqrt(nu)
  return(ri)
}

#'
Initial.points=function(y,X,Z)
{
  X=as.matrix(X)
  Z=as.matrix(Z)
  x2=X[,-1]
  ystar=log(y)-log(1-y)
  lmbeta=suppressWarnings(robustbase::lmrob(ystar~x2))
  betaini=as.numeric(lmbeta$coefficients)
  muini=exp(lmbeta$fitted.values-Rmpfr::log1pexp(lmbeta$fitted.values))
  sigma2ini=stats::sigma(lmbeta)^2*muini^2*(1-muini)^2
  phiini=as.numeric(muini*(1-muini)/sigma2ini)
  if(dim(Z)[2]>1)
  {
    phistar=log(phiini)
    gammaini=solve(t(Z)%*%Z)%*%t(Z)%*%phistar
  }else{
    gammaini=log(mean(phiini,na.rm=T))
  }
  gammaini=as.numeric(c(gammaini))
  Est.param.Rbst=c(betaini,gammaini)
  return(Est.param.Rbst)
}

#'
star.obs=function(p.valor)
{
  obs=NULL
  if(p.valor<0.001){
    obs="***"
  } else if(p.valor<0.01){
    obs="**"
  } else if(p.valor<0.05)
  {
    obs="*"
  } else if(p.valor<0.1)
  {
    obs="."
  } else {
    obs=" "
  }
  return(obs)
}

#'
SQV=function(zq,n,p){
  if(!is.null(zq) && dim(zq)[1]!=1 ){
    return(sqrt(apply((diff(zq,1,1))^2,1,sum))/(sqrt(n)*p))
  }else{
    return(c(0))
  }
}

#'
hatvalues=function(object)
{
  UseMethod("hatvalues")
}
