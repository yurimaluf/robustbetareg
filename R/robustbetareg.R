#' Robust Beta Regression
#'
#' Fit robust beta regression models for rates and proportions via LSMLE, LMDPDE,
#'  SMLE and MDPDE. Both mean and precision of the response variable are modeled
#'  through parametric functions.
#'
#' @name robustbetareg
#'
#' @param formula symbolic description of the model. See Details for further
#'   information.
#' @param data dataset to be used.
#' @param alpha numeric in \eqn{[0,1)} indicating the value of the tuning constant
#'   alpha. \code{alpha = 0} leads to the maximum likelihood estimator.
#'   Robust procedures require \code{alpha} greater than zero.
#'   If this argument is suppressed, the tuning constant will be selected
#'   automatically through the data-driven algorithm proposed by Ribeiro and
#'   Ferrari (2022).
#' @param type character specifying the type of robust estimator to be used in the
#'   estimation process. Supported estimators are "\code{LSMLE}" ,
#'   "\code{LMDPDE}",  "\code{SMLE}", and "\code{MDPDE}"; for details, see Maluf
#'    et al. (2022). The "\code{LSMLE}" is the default.
#' @param link an optional character that specifies the link function of the
#'   mean submodel (mu). The "\code{logit}", "\code{probit}", "\code{cloglog}",
#'   "\code{cauchit}", "\code{loglog}" functions are supported. The \code{logit}
#'   function is the default.
#' @param link.phi an optional character that specifies the link function of the
#'   precision submodel (phi). The "\code{identity}", "\code{log}", "\code{sqrt}"
#'   functions are supported. The default is \code{log} unless formula is of type
#'   \code{y ~ x} where the default is "\code{identity}".
#' @param control a list of control arguments specified via
#'   \code{\link{robustbetareg.control}}.
#' @param model logical. If \code{TRUE} the corresponding components of the fit
#'   (model frame, response, model matrix) are returned.
#' @param y,x,z \code{y} must be a numeric response vector (with values in
#'   \eqn{(0,1)}), \code{x} must be a numeric regressor matrix for the mean
#'    submodel, and \code{z} must be a numeric regressor matrix for the precision
#'    submodel.
#' @param ... argument to be passed to \code{\link{robustbetareg.control}}.
#'
#' @details Beta regression models are employed to model continuous response
#'    variables in the unit interval, like rates and proportions. The maximum
#'    likelihood-based inference suffers from
#'    the lack of robustness in the presence of outliers. Based on
#'    the density power divergence, Ghosh (2019) proposed the minimum density
#'    power divergence estimator (MDPDE). Ribeiro and Ferrari (2022) proposed an
#'    estimator based on the maximization of a reparameterized Lq-likelihood;
#'    it is called SMLE. These estimators require suitable restrictions in the
#'    parameter space. Maluf et al. (2022) proposed robust estimators based on
#'    the MDPDE and the SMLE which have the advantage of overcoming this drawback.
#'    These estimators are called LMDPDE and LSMLE. For details, see the
#'    cited works. The four estimators are implemented in the \code{robustbetareg}
#'    function. They depend on a tuning constant (called \eqn{\alpha}).
#'    When the tuning constant is fixed and equal to 0, all of the estimators
#'    coincide with the maximum likelihood estimator. Ribeiro and Ferrari (2022)
#'    and Maluf et al. (2022) suggest using a data-driven algorithm to select the
#'    optimum value of \eqn{\alpha}. This algorithm is implemented in
#'     \code{robustbetareg} by default when the argument "\code{alpha}" is
#'     suppressed.\cr \cr
#'     The formulation of the model has the same structure as in the usual functions
#'     \code{\link{glm}} and \code{\link{betareg}}. The argument
#'     \code{formula} can comprise of three parts (separated by the symbols
#'     "\eqn{~}" and "\eqn{|}"), namely: observed response variable in the unit
#'     interval, predictor of the mean submodel, with link function \code{link}
#'     and predictor of the precision submodel, with \code{link.phi}
#'     link function. If the model has constant precision, the third part may be
#'     omitted and the link function for phi is "\code{identity}" by default.
#'     The tuning constant \code{alpha} may be treated as fixed or not (chosen
#'     by the data-driven algorithm). If \code{alpha} is fixed, its value
#'     must be specified in the \code{alpha} argument. \cr \cr
#'     Some methods are available for objects of class "\code{robustbetareg}",
#'     see \code{\link{plot.robustbetareg}}, \code{\link{summary.robustbetareg}},
#'     \code{\link{coef.robustbetareg}}, and \code{\link{residuals.robustbetareg}},
#'      for details and other methods.
#'
#'
#' @author Yuri S. Maluf (\email{yurimaluf@@gmail.com}), Francisco F. Queiroz (\email{ffelipeq@@outlook.com}) and Silvia L. P. Ferrari.
#' @references Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust
#'    beta regression through the logit transformation. \emph{arXiv}:2209.11315.\cr \cr
#'    Ribeiro, T.K.A. and Ferrari, S.L.P.  (2022). Robust estimation in beta regression
#'    via maximum Lq-likelihood. \emph{Statistical Papers}. DOI: 10.1007/s00362-022-01320-0. \cr \cr
#'    Ghosh, A. (2019). Robust inference under the beta regression model with
#'    application to health care studies. \emph{Statistical Methods in Medical
#'    Research}, 28:271-888.\cr \cr
#'    Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). On beta regression residuals. \emph{Journal of Applied Statistics}, 35:407–419.
#' @seealso \code{\link{robustbetareg.control}}, \code{\link{summary.robustbetareg}}, \code{\link{residuals.robustbetareg}}
#'
#' @return \code{robustbetareg} returns an object of class "\code{robustbetareg}" with a list of the following components:\tabular{ll}{
#'    \code{coefficients} \tab a list with the "\code{mean}" and "\code{precision}"
#'    coefficients. \cr
#'    \tab \cr
#'    \code{vcov} \tab covariance matrix. \cr
#'    \tab \cr
#'    \code{converged} \tab  logical indicating successful convergence of the
#'       iterative process. \cr
#'    \tab \cr
#'    \code{fitted.values} \tab a vector with the fitted values of the mean submodel. \cr
#'    \tab \cr
#'    \code{start} \tab a vector with the starting values used in the iterative process. \cr
#'    \tab \cr
#'    \code{weights} \tab the weights of each observation in the estimation process. \cr
#'    \tab \cr
#'    \code{Tuning} \tab value of the tuning constant (automatically chosen or fixed) used
#'       in the estimation process. \cr
#'    \tab \cr
#'    \code{residuals} \tab a vector of standardized weighted residual 2 (see Espinheira et al. (2008)). \cr
#'    \tab \cr
#'    \code{n} \tab number of observations. \cr
#'    \tab \cr
#'    \code{link} \tab link function used in the mean submodel. \cr
#'    \tab \cr
#'    \code{link.phi} \tab link function used in the precision submodel. \cr
#'    \tab \cr
#'    \code{Optimal.Tuning} \tab logical indicating whether the data-driven algorithm
#'       was used. \cr
#'    \tab \cr
#'    \code{pseudo.r.squared} \tab pseudo R-squared value. \cr
#'    \tab \cr
#'    \code{control} \tab the control arguments passed to the data-driven algorithm and
#'      \code{optim} call. \cr
#'    \tab \cr
#'    \code{std.error} \tab the standard errors. \cr
#'    \tab \cr
#'    \code{method} \tab type of estimator used. \cr
#'    \tab \cr
#'    \code{call} \tab the original function call. \cr
#'    \tab \cr
#'    \code{formula} \tab the formula used. \cr
#'    \tab \cr
#'    \code{model} \tab the full model frame. \cr
#'    \tab \cr
#'    \code{terms} \tab a list with elements "\code{mean}", "\code{precision}" and "\code{full}"
#'        containing the term objects for the respective models.  \cr
#'    \tab \cr
#'    \code{y} \tab the response variable. \cr
#'    \tab \cr
#'    \code{data} \tab the dataset used. \cr
#' }
#'
#' @examples
#' #### Risk Manager Cost data
#' #data("RiskManagerCost")
#' #
#' ## MLE fit (fixed alpha equal to zero)
#' #fit_MLE <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
#' #                         data = RiskManagerCost, type = "LMDPDE", alpha = 0,
#' #                         link.phi = "log")
#' #summary(fit_MLE)
#' #
#' ## MDPDE with alpha = 0.04
#' #fit_MDPDE <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
#' #                           data = RiskManagerCost, type = "MDPDE",
#' #                           alpha = 0.04, link.phi = "log")
#' #summary(fit_MDPDE)
#' #
#' ## Choosing alpha via data-driven algorithm
#' #fit_MDPDE2 <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
#' #                            data = RiskManagerCost, type = "MDPDE",
#' #                            link.phi = "log")
#' #summary(fit_MDPDE2)
#' #
#' ## Similar result for the LMDPDE fit:
#' #fit_LMDPDE2 <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
#' #                             data = RiskManagerCost, type = "LMDPDE",
#' #                             link.phi = "log")
#' #summary(fit_LMDPDE2)
#'
#' \dontrun{
#' # Diagnostic plots
#' plot(fit_LMDPDE2)
#'
#' #### HIC data
#' data("HIC")
#'
#' # MLE (fixed alpha equal to zero)
#' fit_MLE <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita |
#'                          GDP_percapita, data = HIC, type = "LMDPDE",
#'                          alpha = 0)
#' summary(fit_MLE)
#'
#' # SMLE and MDPDE with alpha selected via data-driven algorithm
#' fit_SMLE <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita |
#'                           GDP_percapita, data = HIC, type = "SMLE")
#' summary(fit_SMLE)
#' fit_MDPDE <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita |
#'                            GDP_percapita, data = HIC, type = "MDPDE")
#' summary(fit_MDPDE)
#' # SMLE and MDPDE return MLE because of the lack of stability
#'
#' # LSMLE and LMDPDE with alpha selected via data-driven algorithm
#' fit_LSMLE <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita |
#'                            GDP_percapita, data = HIC, type = "LSMLE")
#' summary(fit_LSMLE)
#' fit_LMDPDE <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita |
#'                             GDP_percapita, data = HIC, type = "LMDPDE")
#' summary(fit_LMDPDE)
#' # LSMLE and LMDPDE return robust estimates with alpha = 0.06
#'
#'
#' # Plotting the weights against the residuals - LSMLE fit.
#' plot(fit_LSMLE$residuals, fit_LSMLE$weights, pch = "+", xlab = "Residuals",
#'      ylab = "Weights")
#' #identify(fit_LSMLE$residuals, fit_LSMLE$weights) # see observation #1
#'
#' # Excluding outlier observation.
#' fit_LSMLEwo1 <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita |
#'                               GDP_percapita, data = HIC[-1,], type = "LSMLE")
#' summary(fit_LSMLEwo1)
#'
#' # Normal probability plot with simulated envelope
#' plotenvelope(fit_LSMLE)}
#' @import betareg
#' @importFrom stats as.formula model.frame model.response model.matrix terms
#'        delete.response optim qlogis cor var dbeta
#' @rdname robustbetareg
#' @export
robustbetareg = function(formula, data, alpha, type = c("LSMLE","LMDPDE","SMLE","MDPDE"),link = c("logit", "probit", "cloglog", "cauchit", "loglog"),
                         link.phi = NULL,control = robustbetareg.control(...), model = TRUE,...)
{
  cl = match.call()
  type = match.arg(type)
  ocontrol=control
  if(missing(data)){data=environment(formula)}
  if(!missing(alpha)){control$alpha.optimal=FALSE}
  if(missing(alpha)){alpha=NULL}
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf$drop.unused.levels = TRUE
  formula = Formula::as.Formula(formula)
  oformula = as.formula(formula)


  if(length(formula)[2L] < 2L) {
    formula = Formula::as.Formula(formula(formula), ~1)
    simple_formula = TRUE
  }else {
    if(length(formula)[2L] > 2L) {
      formula = Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf1=model.frame(formula,data=data)
  y=model.response(mf1)
  x=model.matrix(formula,data=mf1,rhs = 1L)
  z=model.matrix(formula,data=mf1,rhs = 2L)
  if(simple_formula){colnames(z)[1]="(Phi)"}

  mt = terms(formula, data = data)
  mtX = terms(formula, data = data, rhs = 1L)
  mtZ = delete.response(terms(formula, data = data, rhs = 2L))

  #Error Model Treatment
  if(length(y) < 1){stop("empty model")}
  if(!(min(y) > 0 & max(y) < 1)){stop("invalid dependent variable, all observations must be in (0, 1)")}
  if(!is.null(control$start) & (ncol(x)+ncol(z))!=length(control$start) ){stop("Invalid initial starting point")}
  if(!is.null(alpha)){if(alpha < 0 || alpha > 1){stop("invalid tuning constant, the value must be in [0, 1)")}}
  if(!is.null(link.phi))
  {
    if(link.phi=="identity" & !simple_formula){link.phi="log";warning("Non suitable precision link function, log link used instead")}
  }else{
    link.phi = if(simple_formula){"identity"}
    else "log"
  }
  link = match.arg(link)
  linkobj = set.link(link.mu = link, link.phi = link.phi)

  if(is.null(control$start))
  {
    est.mle=suppressWarnings(betareg(oformula,data,link=link,link.phi = link.phi))
    control$start=c(est.mle$coefficients$mean,est.mle$coefficients$precision)
  }
  if(simple_formula){#Transformation phi -> gamma
    control$start[length(control$start)]=linkobj$linkfun.phi$linkfun(control$start[length(control$start)])
  }

  if(type=="LMDPDE")
  {
    result=LMDPDE.fit(y,x,z,alpha=alpha,link=link,link.phi=link.phi,control=control)
  }
  if(type=="LSMLE")
  {
    result=LSMLE.fit(y,x,z,alpha=alpha,link=link,link.phi=link.phi,control=control)
  }
  if(type=="SMLE")
  {
    result=SMLE.fit(y,x,z,alpha=alpha,link=link,link.phi=link.phi,control=control)
  }
  if(type=="MDPDE")
  {
    result=MDPDE.fit(y,x,z,alpha=alpha,link=link,link.phi=link.phi,control=control)
  }
  result$y=y
  if(simple_formula){
    precision.name=names(result$coefficients$precision)
    result$coefficients$precision=linkobj$linkfun.phi$inv.link(z%*%result$coefficients$precision)[1]
    names(result$coefficients$precision)=precision.name
  }
  if(model){result$model=list(mean = x, precision = z)}
  result$terms=list(mean = mtX, precision = mtZ, full = mt)
  result$call = cl
  result$data = mf1
  result$formula=as.formula(formula)
  gc()
  return(result)
}


#' @rdname robustbetareg
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
  result$weights=dEGB(y_star,mu_hat,phi_hat)^(alpha)#Weights
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
    result$message=c(str1,str2)
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


#' @keywords internal
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

  #Check robustbetareg.control starting points
  if(is.null(control$start)){
    #mle starting point attempt
    est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    if(is.null(est.log.lik))
    {
      est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi, control=betareg.control(start=Initial.points(y,x,z)))),error=function(e) NULL)
    }
    #checking mle starting point
    if(!is.null(est.log.lik)){
      control$start=do.call("c",est.log.lik$coefficients)
      names(control$start)=c(colnames(x),colnames(z))
    }else{
      control$start=Initial.points(y,x,z)
      names(control$start)=c(colnames(x),colnames(z))
    }
  }else{
    names(control$start)=c(colnames(x),colnames(z))
  }
  p=length(control$start)
  control.temp=control

  for(k in 1:M1)#First M1 attempts of tuning parameters
  {
    LMDPDE.par=tryCatch(LMDPDE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control.temp),error=function(e) {LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(LMDPDE.par$converged)
    {
      control.temp$start=c(LMDPDE.par$coefficients$mean,LMDPDE.par$coefficients$precision)
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
  sqv=as.numeric(SQV(zq.t,n,p))
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    LMDPDE.par.star<-tryCatch(LMDPDE.list[[1]],error=function(e){LMDPDE.par.star<-LMDPDE.par;LMDPDE.par.star$message<-"The function cannot be evaluated on initial parameters";return(LMDPDE.par.star)})
    LMDPDE.par.star$sqv=sqv
    LMDPDE.par.star$Optimal.Tuning=TRUE
    rm(LMDPDE.list)
    return(LMDPDE.par.star)
  }
  if(unstable)#Lack of stability
  {
    LMDPDE.par.star<-tryCatch(LMDPDE.list[[1]],error=function(e){LMDPDE.par.star<-LMDPDE.par;return(LMDPDE.par.star)})
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
    LMDPDE.par=tryCatch(LMDPDE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control.temp),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(LMDPDE.par$converged)
    {
      control.temp$start=c(LMDPDE.par$coefficients$mean,LMDPDE.par$coefficients$precision)
    }
    if(any(is.na(do.call("c",LMDPDE.par$coefficients)/do.call("c",LMDPDE.par$std.error))) || is.null(do.call("c",LMDPDE.par$std.error)))
    {
      unstable=T
      break
    }
    LMDPDE.list[[k]]=LMDPDE.par
    zq.t=unname(rbind(zq.t,do.call("c",LMDPDE.par$coefficients)/do.call("c",LMDPDE.par$std.error)))
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


#' @keywords internal
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
    D_q = sum(log(dEGB(y_star,mu_hat,phi_hat)))
  }else{
    a0 = mu_hat*phi_hat
    b0 = (1-mu_hat)*phi_hat
    a_alpha = a0*(1+alpha)
    b_alpha = b0*(1+alpha)
    E_alpha = exp(lgamma(a_alpha)+lgamma(b_alpha)-lgamma(a_alpha+b_alpha)-(1+alpha)*(lgamma(a0)+lgamma(b0)-lgamma(a0+b0)))
    D_q = sum((1+alpha)/(alpha)*dEGB(y_star,mu_hat,phi_hat)^(alpha)-E_alpha)
  }
  return(D_q)
}

#' @keywords internal
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

  Walpha = diag(Phi_Tb*dEGB(y_star,mu_hat,phi_hat)^(alpha))
  Calpha = diag(Phi_Tb*exp(lgamma(a_alpha)+lgamma(b_alpha)-lgamma(a_alpha+b_alpha)-(1+alpha)*(lgamma(a0)+lgamma(b0)-lgamma(a0+b0))))
  #Calpha = diag(Phi_Tb*exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0)))

  result = (1+alpha)*t(X)%*%Phi_Tb%*%(Walpha%*%Ubeta-Calpha%*%E_Ubeta)
  return(result)
}

#' @keywords internal
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

  Walpha = dEGB(y_star,mu_hat,phi_hat)^(alpha)
  Calpha = exp(lgamma(a_alpha)+lgamma(b_alpha)-lgamma(a_alpha+b_alpha)-(1+alpha)*(lgamma(a0)+lgamma(b0)-lgamma(a0+b0)))
  #Calpha = exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0))

  result = (1+alpha)*t(Z)%*%Tg%*%(Walpha%*%Ugamma-Calpha%*%E_Ugamma)
  return(result)
}

#' @keywords internal
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


#' @keywords internal
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

  inv.Lambda=tryCatch(solve(Lambda),error=function(e) {e})
  if(!BBmisc::is.error(inv.Lambda)){
    V=n*inv.Lambda%*%Sigma%*%t(inv.Lambda)
  }else{
    inv.Lambda=MASS::ginv(Lambda)
    V=n*inv.Lambda%*%Sigma%*%t(inv.Lambda)
  }

  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma

  result$Cov=V/n
  result$Std.Error=suppressWarnings(c(sqrt(diag(V/n))))

  return(result)
}

#' @rdname robustbetareg
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
  result$weights=(dEGB(y_star,mu_hat,phi_hat/q))^(alpha)#Weights
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
    result$message=c(str1,str2)
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


#' @keywords internal
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

  #Check robustbetareg.control starting points
  if(is.null(control$start)){
    #mle starting point attempt
    est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    if(is.null(est.log.lik))
    {
      est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi, control=betareg.control(start=Initial.points(y,x,z)))),error=function(e) NULL)
    }
    #checking mle starting point
    if(!is.null(est.log.lik)){
      control$start=do.call("c",est.log.lik$coefficients)
      names(control$start)=c(colnames(x),colnames(z))
    }else{
      control$start=Initial.points(y,x,z)
      names(control$start)=c(colnames(x),colnames(z))
    }
  }else{
    names(control$start)=c(colnames(x),colnames(z))
  }
  p=length(control$start)
  control.temp=control

  for(k in 1:M1)#First M1 attempts of tuning parameters
  {
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control.temp),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(LSMLE.par$converged)
    {
      control.temp$start=c(LSMLE.par$coefficients$mean,LSMLE.par$coefficients$precision)
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
  sqv=as.numeric(SQV(zq.t,n,p))
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    LSMLE.par.star<-tryCatch(LSMLE.list[[1]],error=function(e){LSMLE.par.star<-LSMLE.par;LSMLE.par.star$message<-"The function cannot be evaluated on initial parameters";return(LSMLE.par.star)})
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  if(unstable)#Lack of stability
  {
    LSMLE.par.star<-tryCatch(LSMLE.list[[1]],error=function(e){LSMLE.par.star<-LSMLE.par;return(LSMLE.par.star)})
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    LSMLE.par.star$message="Lack of stability"
    return(LSMLE.par.star)
  }
  if(alpha.ind<8){#Which Tuning satisfy the condition of stability
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
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control.temp),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(LSMLE.par$converged)
    {
      control.temp$start=c(LSMLE.par$coefficients$mean,LSMLE.par$coefficients$precision)
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
    k=suppressWarnings(max(1,min(which(zoo::rollapply(sqv<L,M,sum)==M)))+M+1)
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


#' @keywords internal
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
  f_q_star = dEGB(y_star,mu_hat,phi_q)
  if(alpha==0){
    L_q=sum(log(f_q_star))
  }else{
    L_q=sum((f_q_star^(alpha)-1)/alpha)
  }
  return(L_q)
}

#' @keywords internal
Psi_Beta_LSMLE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  phi_q=phi_hat/q

  aq=mu_hat*phi_q
  bq=(1-mu_hat)*phi_q

  y_star=log(y)-log(1-y)
  mu_star=digamma(aq)-digamma(bq)
  Tb = (link.model$link.mu$d.linkfun(mu_hat))^(-1)
  f_q_star = (dEGB(y_star,mu_hat,phi_q))^(alpha)

  result=t(X)%*%diag(phi_q*Tb*f_q_star)%*%(y_star-mu_star)
  return(result)
}

#' @keywords internal
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

  f_q_star = (dEGB(y_star,mu_hat,phi_q))^(alpha)
  eta = mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger)
  result=t(Z)%*%diag(Tg%*%f_q_star/q)%*%eta
  return(result)
}

#' @keywords internal
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

#' @keywords internal
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

  #V=n*inv.Lambda%*%Sigma%*%t(inv.Lambda)
  #V=n*MASS::ginv(Lambda)%*%Sigma%*%t(MASS::ginv(Lambda))
  inv.Lambda=tryCatch(solve(Lambda),error=function(e) {e})
  if(!BBmisc::is.error(inv.Lambda)){
    V=n*inv.Lambda%*%Sigma%*%t(inv.Lambda)
  }else{
    inv.Lambda=MASS::ginv(Lambda)
    V=n*inv.Lambda%*%Sigma%*%t(inv.Lambda)
  }

  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma
  result$Cov=V/n

  result$K=t(X)%*%Tb%*%Q.inv%*%(Phi^2)%*%C.k0%*%Lambda_mu_mu%*%Tb%*%X
  result$Std.Error=suppressWarnings(c(sqrt(diag(V/n))))

  return(result)
}

#' @rdname robustbetareg
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

  critical.area=NULL
  theta=tryCatch(optim(par=start_theta,fn=D_q,gr=Psi_MDPDE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1)),error=function(e){
    theta$msg<-e$message;check<<-F;return(theta)},
    warning=function(w){
      critical.area<<-w$message;
      theta<-suppressWarnings(optim(par=start_theta,fn=D_q,gr=Psi_MDPDE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1)));return(theta)})

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
  str1=str2=str3=NULL
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
  if(length(critical.area)!=0)
  {
    str3="In some optimization interactions, one or more points reached the critical area."
  }
  if(!is.null(str1)||!is.null(str2))
  {
    result$message=c(str1,str2)
  }
  if(!is.null(str3)){
    result$warning=c(str3)
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

#' @keywords internal
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

  #Check robustbetareg.control starting points
  if(is.null(control$start)){
    #mle starting point attempt
    est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    if(is.null(est.log.lik))
    {
      est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi, control=betareg.control(start=Initial.points(y,x,z)))),error=function(e) NULL)
    }
    #checking mle starting point
    if(!is.null(est.log.lik)){
      control$start=do.call("c",est.log.lik$coefficients)
      names(control$start)=c(colnames(x),colnames(z))
    }else{
      control$start=Initial.points(y,x,z)
      names(control$start)=c(colnames(x),colnames(z))
    }
  }else{
    names(control$start)=c(colnames(x),colnames(z))
  }
  p=length(control$start)

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
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    MDPDE.par.star<-tryCatch(MDPDE.list[[1]],error=function(e){MDPDE.par.star<-MDPDE.par;MDPDE.par.star$message<-"The function cannot be evaluated on initial parameters";return(MDPDE.par.star)})
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
    return(MDPDE.par.star)
  }
  if(unstable)#Lack of stability
  {
    MDPDE.par.star<-MDPDE.list[[1]]
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
    MDPDE.par.star$message="Lack of stability"
    return(MDPDE.par.star)
  }
  if(alpha.ind<8){#Which Tuning satisfy the condition os stability
    MDPDE.par.star<-MDPDE.list[[alpha.ind+1]]
    MDPDE.par.star$sqv=sqv
    MDPDE.par.star$Optimal.Tuning=TRUE
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

#' @keywords internal
D_q=function(theta,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  k=ncol(X)
  m=ncol(Z)
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)
  Beta=theta[1:k]
  Gamma=theta[(k+1):(k+m)]

  mu_hat = link.model$linkfun.mu$inv.link(X%*%Beta)
  phi_hat = link.model$linkfun.phi$inv.link(Z%*%Gamma)

  if(any(mu_hat==0 || mu_hat==1)){
    mu_hat=pmax(pmin(mu_hat,1-.Machine$double.eps),.Machine$double.eps)
  }

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

#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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
}

#' @rdname robustbetareg
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
  critical.area=NULL
  theta=tryCatch(optim(par=start_theta,fn=L_q,gr=Psi_SMLE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1)),error=function(e){
    theta$msg<-e$message;check<<-F;return(theta)},
    warning=function(w){
      critical.area<<-w$message;
      theta<-suppressWarnings(optim(par=start_theta,fn=L_q,gr=Psi_SMLE,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control = list(fnscale=-1)));return(theta)})

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
  str1=str2=str3=NULL
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
  if(length(critical.area)!=0)
  {
    str3="In some optimization interactions, one or more points reached the critical area."
  }

  if(!is.null(str1)||!is.null(str2))
  {
    result$message=c(str1,str2)
  }
  if(!is.null(str3)){
    result$warning=c(str3)
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

#' @keywords internal
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

  #Check robustbetareg.control starting points
  if(is.null(control$start)){
    #mle starting point attempt
    est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    if(is.null(est.log.lik))
    {
      est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi, control=betareg.control(start=Initial.points(y,x,z)))),error=function(e) NULL)
    }
    #checking mle starting point
    if(!is.null(est.log.lik)){
      control$start=do.call("c",est.log.lik$coefficients)
      names(control$start)=c(colnames(x),colnames(z))
    }else{
      control$start=Initial.points(y,x,z)
      names(control$start)=c(colnames(x),colnames(z))
    }
  }else{
    names(control$start)=c(colnames(x),colnames(z))
  }
  p=length(control$start)

  for(k in 1:M1)#First M1 attempts of tuning parameters
  {
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
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: All M1 sqv beneath L
  {
    SMLE.par.star<-tryCatch(SMLE.list[[1]],error=function(e){SMLE.par.star<-SMLE.par;SMLE.par.star$message<-"The function cannot be evaluated on initial parameters";return(SMLE.par.star)})
    SMLE.par.star$sqv=sqv
    SMLE.par.star$Optimal.Tuning=TRUE
    return(SMLE.par.star)
  }
  if(unstable)#Lack of stability
  {
    SMLE.par.star<-tryCatch(SMLE.list[[1]],error=function(e){SMLE.par.star<-SMLE.par;return(SMLE.par.star)})
    SMLE.par.star$sqv=sqv
    SMLE.par.star$Optimal.Tuning=TRUE
    SMLE.par.star$message="Lack of stability"
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
    k=suppressWarnings(max(1,min(which(zoo::rollapply(sqv<L,M,sum)==M)))+M+1)
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

#' @keywords internal
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
  if(any(mu==0 || mu==1)){
    mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
  }
  a <- mu*phi
  b <- (1 - mu)*phi
  log_likS <- sum(dbeta(y, a, b, log = F)^(alpha))#function to be maximized

  return(log_likS)
}

#' @keywords internal
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

#' @keywords internal
Psi_Gamma_SMLE=function(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi){
  q=1-alpha
  link.model=set.link(link.mu=link_mu,link.phi=link_phi)

  phi_q=q*(phi_hat-2)+2
  mu_q=(q*(mu_hat*phi_hat-1)+1)/phi_q
  y_dagger=log(1-y)
  y_star=log(y)-y_dagger

  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  mu_dagger=digamma((1-mu_hat)*phi_hat)-digamma(phi_hat)

  Tg = (link.model$linkfun.phi$d.linkfun(phi_q))^1

  result=t(Z)%*%diag(mu_q*Tg/q)%*%(mu_q*(y_star-mu_star)+(y_dagger-mu_dagger))
  return(result)
}

#' @keywords internal
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

#' @keywords internal
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
  inv.Jq=tryCatch(solve(Jq), error=function(e) {e})
  if(!BBmisc::is.error(inv.Jq)){
    Vq <- inv.Jq%*%Kq%*%t(inv.Jq)
  }else{
    inv.Jq = MASS::ginv(Jq)
    Vq <- inv.Jq%*%Kq%*%t(inv.Jq)
  }

  result=list()
  result$Lambda=Jq
  result$Sigma=Kq
  result$Cov=Vq
  result$Std.Error=suppressWarnings(t(sqrt(diag(Vq))))

  return(result)
}



#' Auxiliary for Controlling robustbetareg Fitting
#'
#' Several parameters that control fitting of robust beta regression models using
#'  \code{\link[=robustbetareg]{robustbetareg.}}
#'
#' @name robustbetareg.control
#'
#' @param start an optional vector with starting values for the parameter estimates.
#' @param alpha.optimal a logical value. If \code{TRUE} the tuning constant will
#'      be select via the data-driven algorithm.
#' @param tolerance numeric tolerance for convergence.
#' @param maxit argument passed to \code{\link{optim}}.
#' @param L numeric specifying the threshold for the data-driven algorithm.
#'      Default is \code{L = 0.02}.
#' @param M integer specifying the number of grid spacing for the data-driven
#'      algorithm. Default is \code{M = 3}.
#' @param ... currently not used.
#'
#' @details The \code{robustbetareg.control} controls the fitting process of
#'     the robust estimation in beta regression via the LMDPDE, LSMLE, MDPDE, and
#'     SMLE. The arguments \code{L} and \code{M} are passed to the data-driven
#'     algorithm for selecting the optimum alpha value; details can be found in
#'     Ribeiro and Ferrari (2022). Starting values for the parameters associated
#'     to the mean and precision submodels may be supplied via \code{start}. \cr
#'     We do not recommend changing the arguments passed to the data-driven algorithm.
#' @author Yuri S. Maluf (\email{yurimaluf@@gmail.com}),
#' Francisco F. Queiroz (\email{ffelipeq@@outlook.com}) and Silvia L. P. Ferrari.
#'
#' @references Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust
#'    beta regression through the logit transformation. \emph{arXiv}:2209.11315.\cr \cr
#'    Ribeiro, K.A.T. and Ferrari, S.L.P.  (2022). Robust estimation in beta regression
#'    via maximum Lq-likelihood. \emph{Statistical Papers}. DOI: 10.1007/s00362-022-01320-0. \cr \cr
#'
#' @return A list with components named as the arguments.
#' @seealso \code{\link{robustbetareg}}
#' @examples
#' \dontrun{
#' data("RiskManagerCost")
#'
#' # Using a robust start value for the parameters associated with the
#' # mean submodel
#' # using the robustbase package
#' # robust regression to obtain a starting value for beta
#' fit_lm_Rob <- robustbase::lmrob(FIRMCOST ~ SIZELOG + INDCOST,
#'                                 data = RiskManagerCost)
#' initials_beta_rob <-  as.numeric(coef(fit_lm_Rob))
#' etarob <- model.matrix(fit_lm_Rob)%*%initials_beta_rob
#' muhat_Rob <- set.link(link.mu = "logit",
#'                       link.phi = "log")$linkfun.mu$inv.link(etarob)
#' T_1_Rob <- 1/set.link(link.mu = "logit",
#'                       link.phi = "log")$linkfun.mu$d.linkfun(muhat_Rob)
#' #estimate of variance of y based on robustbase package
#' sigma2hat_Rob <- ((fit_lm_Rob$scale^2)*(T_1_Rob^2))
#' #phi estimate from robust method
#' phihat_Rob <- mean((muhat_Rob*(1-muhat_Rob))/sigma2hat_Rob)
#' gama1hat_rob <- log(phihat_Rob)
#' #gamma estimates from robustbase package
#' initials_gama_rob <-  as.numeric(gama1hat_rob)
#' #robust starting values for beta and gamma
#' thetaStart <- c(initials_beta_rob, initials_gama_rob)
#'
#' fit_LSMLE <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
#'                            data = RiskManagerCost,
#'                            type = "LSMLE", link.phi = "log",
#'                            control = robustbetareg.control(start = thetaStart))
#' }
#' @export
robustbetareg.control=function(start = NULL, alpha.optimal = TRUE, tolerance = 1e-3,
                               maxit = 5000, L = 0.02, M = 3, ...)
{
  UseMethod("robustbetareg.control")
}

#' @export
robustbetareg.control.default=function(start=NULL,alpha.optimal=TRUE,tolerance=1e-3,maxit=5000,L=0.02,M=3,...)
{
  result <- list(start=start,alpha.optimal=alpha.optimal,tolerance = tolerance, maxit = maxit, L=L, M=M)
  result <- c(result, list(...))
  return(result)
}


#' The Exponential Generalized Beta of the Second Type Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for exponential generalized beta of the second type distribution.
#'
#' @name EGB
#' @param y_star,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#' to be the number required.
#' @param mu mu parameter.
#' @param phi phi parameter.
#' @param log logical; if \code{TRUE}, probabilities p are given as log(p). Default is FALSE.
#'
#' @details The EGB distribution with parameters \code{mu = }\eqn{\mu} and
#'     \code{phi = }\eqn{\phi} has density  \deqn{f(y^\star;\mu,\phi)=
#'     B^{-1}(\mu\phi,(1-\mu)\phi) \exp\{-y^\star(1-\mu)\phi\}/ (1+\exp\{-y^\star\})^{\phi},}
#'     with \eqn{\mu\in(0,1),\phi>0} and \eqn{y^\star \in (-\infty, \infty)}. For this
#'     distribution, \eqn{E(y^\star)=\psi(\mu\phi)-\psi((1-\mu)\phi)} and
#'     \eqn{Var(y^\star)=\psi'(\mu\phi)+\psi'((1-\mu)\phi)}, where \eqn{\psi}
#'     is the digamma function. See Kerman and McDonald (2015) for additional
#'     details. If \eqn{y \sim beta(\mu, \phi)}, with \eqn{\mu} and
#'     \eqn{\phi} representing the mean and precision of \eqn{y}, then
#'     \eqn{y^\star = \log(y/(1-y)) \sim EGB(\mu, \phi)} with the density
#'     given above.
#' @author Yuri S. Maluf (\email{yurimaluf@@gmail.com}),
#' Francisco F. Queiroz (\email{ffelipeq@@outlook.com}) and Silvia L. P. Ferrari.
#'
#' @references Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust
#'     beta regression through the logit transformation. \emph{arXiv}:2209.11315.\cr \cr
#'     Kerman, S. and McDonald, J.B. (2015). Skewness-kurtosis bounds for EGB1, EGB2,
#'     and special cases. \emph{Communications in Statistics - Theory and Methods},
#'     44:3857-3864.
#'
#' @examples
#' dEGB(0.2, mu = 0.3, phi = 1, )
#' mu = 0.2; phi = 2;
#' set.seed(1)
#' EGBsample = rEGB(1000, mu, phi)
#' hist(EGBsample, prob = TRUE, breaks = 15, main = "", las = 1, ylim = c(0, 0.2),
#'      xlim = c(-20, 10))
#' curve(dEGB(x, mu, phi), from = -20, to = 8, add = TRUE, col = "red")
#' rug(EGBsample)
#'
#'
#' # Showing the P(Y* < -5) = 0.17, where Y* ~ EGB(0.2, 2).
#' x = seq(-20, 10,0.01)
#' y = dEGB(x, mu, phi)
#' plot(x, y, type = "l", lwd = 2, las = 1)
#' x1 = seq(-20, -5, 0.01)
#' y1 = dEGB(x1, mu, phi)
#' polygon(c(x1, -5, -5), c(y1, 0, 0), col = "lightblue")
#'
#'
#' plot(x, pEGB(x, mu, phi), type = "l", las = 1, lwd = 2,
#'      ylab = expression(P("Y*"<y)), xlab = "y")
#' p = pEGB(0, mu, phi)
#' q = qEGB(p, mu, phi)
#' points(q, p, pch = 16, col = 2, cex = 1.5)
#' text(2, 0.83, paste("(", 0, ",", round(p, 2), ")"), font = 2,
#'      cex = 0.8, col = "red")
#'
#' @return \code{dEGB} gives the density, \code{pEGB} gives the distribution function,
#'     \code{qEGB} gives the quantile function, and \code{rEGB} generates random
#'      variables.
#' @importFrom stats pbeta qbeta rbeta
#'@export
dEGB=function(y_star, mu, phi, log = FALSE)
{
  #options(warn = 2) #Convert warnings into errors
  if (any((-abs(2*mu-1)+1)<=0)){
    return(warning("'mu' parameter must be within unit interval"))
  }
  if (any(phi<=0)){
    return(warning("'phi' parameter must be a positive value"))
  }

  a0=mu*phi
  b0=(1-mu)*phi
  if(log==F){
    k=tryCatch(exp(-(lgamma(a0)+lgamma(b0)-lgamma(a0+b0)+b0*y_star+phi*Rmpfr::log1pexp(-y_star))),error=function(e){stop("Error")})
  }
  if(log==T){
    k=tryCatch(-(lgamma(a0)+lgamma(b0)-lgamma(a0+b0)+b0*y_star+phi*Rmpfr::log1pexp(-y_star)),error=function(e){stop("Error")})
  }
  return(k)
}

#' @rdname EGB
#' @export
pEGB=function(q, mu, phi)
{
  if (any((-abs(2*mu-1)+1)<=0)){
    return(warning("'mu' parameter must be within unit interval"))
  }
  if (any(phi<=0)){
    return(warning("'phi' parameter must be a positive value"))
  }
  return(pbeta((1+exp(-q))^(-1),mu*phi,(1-mu)*phi))
}


#' @rdname EGB
#' @export
qEGB=function(p, mu, phi)
{
  if (any((-abs(2*mu-1)+1)<=0)){
    return(warning("'mu' parameter must be within unit interval"))
  }
  if (any(phi<=0)){
    return(warning("'phi' parameter must be a positive value"))
  }
  q=qbeta(p,mu*phi,(1-mu)*phi)
  return(-log((1/q)-1))
}

#' @rdname EGB
#' @export
rEGB=function(n, mu, phi)
{
  if (any((-abs(2*mu-1)+1)<=0)){
    return(warning("'mu' parameter must be within unit interval"))
  }
  if (any(phi<=0)){
    return(warning("'phi' parameter must be a positive value"))
  }
  h=rbeta(n,mu*phi,(1-mu)*phi)
  return(log(h)-log(1-h))
}

#' Normal Probability Plots of Residuals with Simulated Envelope for robustbetareg Objects
#'
#' \code{plotenvelope} is used to display normal probability plots of residuals
#' with simulated envelope for the robust beta regression. Currently, eight types of
#' residuals are supported: sweighted2, pearson, weighted, sweighted,
#' sweighted.gamma, sweighted2.gamma, combined, and combined.projection residuals.
#'
#' @name plotenvelope
#'
#' @param object fitted model object of class \code{robustbetareg}.
#' @param type character indicating the type of residuals to be used, see
#'      \code{\link{residuals.robustbetareg}}. Default is \code{type = "sweighted2"}.
#' @param conf numeric specifying the confidence level of the simulated
#'      envelopes. Default is \code{conf = 0.95}.
#' @param n.sim a positive integer representing the number of iterations
#'      to generate the simulated envelopes. Default is \code{n.sim=100}.
#' @param PrgBar logical. If \code{PrgBar = TRUE} the progress bar will be shown
#'      in the console. Default is \code{PrgBar = TRUE}.
#' @param control a list of control arguments specified via
#'      \code{\link[robustbetareg:robustbetareg.control]{robustbetareg.control}}.
#' @param ... arguments passed to \code{\link{plot}}.
#'
#' @details The \code{plotenvelope} creates normal probability plots with simulated
#'      envelope (see Atkinson (1985) for details). Under the correct model,
#'      approximately 100*conf of the residuals are expected to be inside the
#'      envelope.
#' @return \code{plotenvelope} returns normal probability plot of residuals with simulated
#'      envelope.
#'
#' @author Yuri S. Maluf (\email{yurimaluf@@gmail.com}),
#' Francisco F. Queiroz (\email{ffelipeq@@outlook.com}) and Silvia L. P. Ferrari.
#'
#' @references Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust
#'     beta regression through the logit transformation. \emph{arXiv}:2209.11315.\cr \cr
#'     Atkinson, A.C. (1985) Plots, transformations and regression: an
#'     introduction to graphical methods of diagnostic regression analysis.
#'     \emph{Oxford Science Publications}, Oxford.
#'
#' @seealso \code{\link[robustbetareg:robustbetareg]{robustbetareg}}, \code{\link[robustbetareg:robustbetareg.control]{robustbetareg.control}}, \code{\link[robustbetareg:residuals]{residuals.robustbetareg}}
#'
#' @examples
#' \dontrun{
#' get(data("HIC", package = "robustbetareg"))
#' hic <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita | GDP_percapita,
#' data = HIC, alpha = 0.06)
#' plotenvelope(hic, n.sim = 500)
#'
#' get(data("RiskManagerCost", package = "robustbetareg"))
#' rmc <- robustbetareg(FIRMCOST ~ INDCOST + SIZELOG | INDCOST + SIZELOG, data = RiskManagerCost)
#' plotenvelope(rmc, conf = 0.90)}
#'
#' @export
plotenvelope=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3"),conf=0.95,n.sim=100,PrgBar=T,control=robustbetareg.control(...), ...)
{
  UseMethod("plotenvelope")
}

#' @export
plotenvelope.robustbetareg=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3"),conf=0.95,n.sim=100,PrgBar=T,control=robustbetareg.control(...), ...)
{
  if(missing(control)){control=robustbetareg.control(object)}
  type = match.arg(type)
  ylim.boolean=methods::hasArg(ylim)
  arg=list(...)
  limit=FALSE
  y.sim=ResEnvelop=NULL
  link=object$link
  link.phi=object$link.phi
  y=object$y
  x=as.matrix(object$model$mean)
  z=as.matrix(object$model$precision)
  n=length(object$fitted.values$mu.predict)
  a=object$fitted.values$mu.predict*object$fitted.values$phi.predict
  b=(1-object$fitted.values$mu.predict)*object$fitted.values$phi.predict
  residual=residuals(object,type=type)
  k=1
  if(PrgBar){pb = utils::txtProgressBar(min = 0, max = n.sim, style = 3)}
  while(k<=n.sim & !limit)
  {
    y.sim=pmax(pmin(sapply(seq(1,n,1),function(i) rbeta(1,a[i],b[i])),1-.Machine$double.eps),.Machine$double.eps)
    est.mle=betareg.fit(x,y,z,link=link,link.phi = link.phi)
    start=c(est.mle$coefficients$mean,est.mle$coefficients$precision)
    if(object$method=="LSMLE"){
      robustbetareg.sim=tryCatch(LSMLE.fit(y=y.sim,x=x,z=z,alpha=object$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){robustbetareg.sim$converged<-FALSE; return(robustbetareg.sim)})
    }else if(object$method=="LMDPDE"){
      robustbetareg.sim=tryCatch(LMDPDE.fit(y=y.sim,x=x,z=z,alpha=object$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){robustbetareg.sim$converged<-FALSE; return(robustbetareg.sim)})
    }else if(object$method=="SMLE"){
      robustbetareg.sim=tryCatch(SMLE.fit(y=y.sim,x=x,z=z,alpha=object$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){robustbetareg.sim$converged<-FALSE; return(robustbetareg.sim)})
    }else{
      robustbetareg.sim=tryCatch(MDPDE.fit(y=y.sim,x=x,z=z,alpha=object$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){robustbetareg.sim$converged<-FALSE; return(robustbetareg.sim)})
    }
    if(robustbetareg.sim$converged)
    {
      if(type=="sweighted2")
      {
        ResEnvelop=rbind(ResEnvelop,sort(robustbetareg.sim$residuals,decreasing = F))
      }else{

        robustbetareg.sim$y=y.sim
        robustbetareg.sim$model$mean=object$model$mean
        robustbetareg.sim$model$precision=object$model$precision
        ResEnvelop=rbind(ResEnvelop,sort(residuals(robustbetareg.sim,type=type),decreasing = F))
      }
      k=k+1
      # update progress bar
      if(PrgBar){utils::setTxtProgressBar(pb, k)}
    }
    if(k==(2*n.sim)){limit=T}
  }
  Envelope=apply(ResEnvelop,2,stats::quantile,c((1-conf)/2,0.5,1-(1-conf)/2))
  if(!ylim.boolean)
  {
    ylim <- range(Envelope[1,],Envelope[2,],Envelope[3,],residual)
    arg$ylim=ylim
  }
  graphics::par(mar=c(5.0,5.0,4.0,2.0),pch=16, cex=1.0, cex.lab=1.0, cex.axis=1.0, cex.main=1.5)
  #ARG=append(list(y=residual,main="Envelope Plot", xlab="Normal quantiles",ylab="Residuals"),arg)
  ARG=append(list(y=residual,main="", xlab="Normal quantiles",ylab="Residuals"),arg)
  do.call(stats::qqnorm,ARG)
  graphics::par(new=T)
  ARG=utils::modifyList(ARG,list(y=Envelope[1,],axes=F,main = "",xlab="",ylab="",type="l",lty=1,lwd=1.0))
  do.call(stats::qqnorm,ARG)
  graphics::par(new=T)
  ARG=utils::modifyList(ARG,list(y=Envelope[2,],lty=2))
  do.call(stats::qqnorm,ARG)
  graphics::par(new=T)
  ARG=utils::modifyList(ARG,list(y=Envelope[3,],lty=1))
  do.call(stats::qqnorm,ARG)
}


#' Robust Wald-type Tests
#'
#' \code{waldtypetest} provides Wald-type tests for both simple and composite
#' hypothesis for independent but non-homogeneous observations,
#' based on the robust estimators (LSMLE, LMDPDE, SMLE, and MDPDE).
#'
#' @name waldtypetest
#' @param object fitted model object of class \code{robustbetareg} (see \code{\link[robustbetareg:robustbetareg]{robustbetareg}}).
#' @param FUN function representing the null hypothesis to be tested.
#' @param ... further arguments to be passed to the \code{FUN} function.
#'
#' @details The function provides a procedure to test a general hypothesis
#'     \eqn{m(\theta) = \eta_0}, for a fixed \eqn{\eta_0 \in R^d}, against
#'     a two side alternative, through a robust Wald-type test; see Maluf
#'     et al. (2022) for further details. The argument \code{FUN} specifies the
#'     function \eqn{m(\theta) - \eta_0} which defines the null hypothesis to be
#'     considered in the test. For instance, consider a model with
#'     \eqn{\theta = (\beta_1, \beta_2, \beta_3, \gamma_1)^\top} and suppose that we
#'     want to test the null hypothesis \eqn{\beta_2 + \beta_3 = 4}
#'     against a two side alternative. The \code{FUN} argument can be
#'     \code{FUN = function(theta) \{theta[2] + theta[3] - 4\} }. It is also possible
#'     define the function as \code{FUN = function(theta, B) \{theta[2] + theta[3] - B\}},
#'     and pass the \code{B} argument to the \code{waldtypetest} function.
#'     If no function is specified, that is, \code{FUN=NULL}, the \code{waldtypetest}
#'     returns a test in which the null hypothesis is that all the coefficients are
#'     jointly equal to zero.
#'
#' @return \code{waldtypetest} returns an output for the Wald-type test containing
#' the value for the test statistics, degrees-of-freedom and p-value.
#'
#' @references Maluf, Y. S., Ferrari, S. L. P., and Queiroz, F. F. (2022). Robust
#'     beta regression through the logit transformation. \emph{arXiv}:2209.11315.\cr \cr
#'     Basu, A., Ghosh, A., Martin, N., and Pardo, L. (2018). Robust Wald-type
#'     tests for  non-homogeneous observations based on the minimum density
#'     power divergence estimator. \emph{Metrika}, 81:493–522. \cr \cr
#'     Ribeiro, K. A. T. and Ferrari, S. L. P. (2022). Robust estimation in beta
#'     regression via maximum Lq-likelihood. \emph{Statistical Papers}.
#' @author Yuri S. Maluf (\email{yurimaluf@@gmail.com}),
#' Francisco F. Queiroz (\email{ffelipeq@@outlook.com}) and Silvia L. P. Ferrari.
#'
#' @seealso \code{\link[robustbetareg:robustbetareg]{robustbetareg}}
#'
#' @examples
#' \dontrun{
#' set.seed(2022)
#' N <- 40 #Sample Size
#' beta.coef <- c(-1,-2) #Arbitrary Beta Coefficients
#' gamma.coef <- c(5) #Arbitrary Gamma Coefficient
#' X <- cbind(rep(1,N), x <- runif(N))
#' mu <- exp(X%*%beta.coef)/(1+exp(X%*%beta.coef)) #Inverse Logit Link Function
#' phi <- exp(gamma.coef) #Inverse Log Link Function
#' y <- rbeta(N, mu*phi, (1-mu)*phi)
#' y[26] <- rbeta(1,((1 + mu[26])/2)*phi,(1-((1 + mu[26])/2))*phi) #Contaminated data point
#' SimData <- as.data.frame(cbind(y,x))
#' colnames(SimData) <- c("y","x")
#' fit.mle <- robustbetareg(y ~ x | 1, data = SimData, alpha = 0) #Non-Robust Estimator
#' fit.lsmle <- robustbetareg(y ~ x | 1, data = SimData) #Robust Estimator
#' h0 <- function(theta,B){theta[1:2] - B} #Hiphothesis to be tested
#' waldtypetest(fit.mle, h0, B = beta.coef) #Testing beta.1=-1 and beta.2=-2
#' waldtypetest(fit.simdata, h0, B = beta.coef) #Testing beta.1=-1 and beta.2=-2}
#'
#' @importFrom stats pchisq
#'
#' @export
waldtypetest=function(object,FUN,...)
{
  UseMethod("waldtypetest")
}

#' @export
waldtypetest.robustbetareg=function(object,FUN,...)
{
  general=FALSE
  if(missing(FUN)){general=T}
  if(!object$converged){stop(paste("There is no convergence in the model",deparse(substitute(object))))}
  cl = match.call()
  n=length(object$fitted.values$mu.predict)
  V=object$vcov
  b=object$coefficient$mean
  g=object$coefficient$precision
  p=c(b,g)
  k1=length(b)
  k2=length(g)
  k3=length(p)
  if(general)
  {
    result.beta=result.gamma=NULL
    if(ncol(object$model$mean)>1)
    {
      #Beta
      result.beta=list()
      f.b=function(x){x[2:k1]}
      M=numDeriv::jacobian(f.b,p,...)
      m=f.b(p,...)
      r=length(m)
      if(Matrix::rankMatrix(M)[1]!=r)
      {
        stop("The rank matrix is not supported")
      }
      W_alpha=t(m)%*%solve(M%*%V%*%t(M))%*%(m)
      #Beta Register
      result.beta$W.alpha=as.numeric(W_alpha)
      result.beta$df=r
      result.beta$pValue=as.numeric(1-pchisq(W_alpha,df=r))
    }
    if(ncol(object$model$precision)>1)
    {
      result.gamma=list()
      f.g=function(x){x[(k1+2):k3]}
      M=numDeriv::jacobian(f.g,p,...)
      m=f.g(p,...)
      r=length(m)
      if(Matrix::rankMatrix(M)[1]!=r){stop("The rank matrix is not supported")}
      W_alpha=t(m)%*%solve(M%*%V%*%t(M))%*%(m)
      result.gamma$W.alpha=as.numeric(W_alpha)
      result.gamma$df=r
      result.gamma$pValue=as.numeric(1-pchisq(W_alpha,df=r))
    }
    result=list(beta.wald=result.beta,gamma.wald=result.gamma,general=general,msg=paste("Results based on",object$method))
  }else{
    result=list()
    f=FUN
    M=numDeriv::jacobian(f,p,...)
    m=f(p,...)
    r=length(m)
    n=length(object$fitted.values$mu.predict)
    if(Matrix::rankMatrix(M)[1]!=r){stop("The rank matrix is not supported")}
    W_alpha=t(m)%*%solve(M%*%V%*%t(M))%*%(m)
    #Register
    result$W.alpha=as.numeric(W_alpha)
    result$df=r
    result$pValue=as.numeric(1-pchisq(W_alpha,df=r))
    result$general=general
    result$msg=paste("Results based on",object$method)

  }
  result$method=object$method
  class(result)="WaldTest_robustbetareg"
  return(result)
}


#' @export
residuals=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3"),...)
{
  UseMethod("residuals")
}

#' Residuals Method for robustbetareg Objects
#'
#' The function provides several types of residuals for the robust beta regression
#' models: Pearson residuals (raw residuals scaled by square root of variance function)
#' and different kinds of weighted residuals proposed by Espinheira et al. (2008)
#' and Espinheira et al. (2017).
#'
#' @name residuals.robustbetareg
#'
#' @param object fitted model object of class \code{robustbetareg}.
#' @param type character indicating type of residuals to be used.
#' @param ... currently not used.
#'
#' @return \code{residuals} returns a vector with the residuals of the type
#'      specified in the \code{type} argument.
#'
#' @details The definitions of the first four residuals are provided in
#'      Espinheira et al. (2008):  equation (2) for "\code{pearson}",
#'      equation (6) for "\code{weighted}", equation (7) for "\code{sweighted}",
#'      and equation (8) for "\code{sweighted2}". For the last four residuals
#'      the definitions are described in Espinheira et al. (2017): equations (7)
#'      and (10) for the "\code{sweighted.gamma}" and "\code{sweighted2.gamma}",
#'      respectively, equation (9) for "\code{combined}", and equation (11)
#'      for "\code{combined.projection}".
#'
#' @references Maluf, Y. S., Ferrari, S. L. P., and Queiroz, F. F. (2022). Robust
#'     beta regression through the logit transformation. \emph{arXiv}:2209.11315.\cr \cr
#'     Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). On Beta
#'     Regression Residuals. \emph{Journal of Applied Statistics}, 35:407–419.\cr \cr
#'     Espinheira, P.L., Santos, E.G.and Cribari-Neto, F. (2017). On nonlinear
#'     beta regression residuals. \emph{Biometrical Journal}, 59:445-461.\cr \cr
#'
#' @seealso \code{\link[robustbetareg:robustbetareg]{robustbetareg}}
#'
#' @examples
#' \dontrun{
#' get(data("HIC", package = "robustbetareg"))
#' fit.hic <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita | 1, data = HIC, alpha = 0.04)
#' res <- residuals(fit.hic, type = "sweighted2")
#' plot(res)
#' abline(h = 0)}
#'
#' @aliases residuals residuals .robustbetareg
#'
#' @method residuals robustbetareg
#'
#' @export
residuals.robustbetareg=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3"),...)
{
  type = match.arg(type)
  y=object$y
  x=object$model$mean
  z=object$model$precision
  linkobj=set.link(link.mu=object$link,link.phi=object$link.phi)
  mu.predict=object$fitted.values$mu.predict
  phi.predict=object$fitted.values$phi.predict
  switch(type,sweighted2={
    res=sweighted2_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,X=x,linkobj=linkobj)
  },
  sweighted={
    res=res=sweighted_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  pearson={
    res=pearson_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  weighted={
    res=weighted_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  sweighted.gamma={
    res=sweighted.gamma_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  sweighted2.gamma={
    res=sweighted2.gamma_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,Z=z,linkobj=linkobj)
  },
  combined={
    res=combined_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  combined.projection={
    res=combined.projection_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,X=x,Z=z,linkobj=linkobj)
  },
  sweighted3={
    res=sweighted3_res(mu_hat=mu.predict,phi_hat=phi.predict,alpha=object$Tuning,y=y,X=x,linkobj=linkobj)
  },

  stop(gettextf("%s residual not recognised", sQuote(type)),domain = NA))

  return(res)
}


#' Prediction Methods for robustbetareg Objects Class
#'
#' Extract various types of predictions from beta regression models: either on
#' the scale of responses in (0, 1) or the scale of the linear predictor,
#' from \code{robustbetareg} objects.
#'
#'@name predict
#'
#' @param object fitted model object of class "\code{robustbetareg}".
#' @param newdata optional, a data frame with new predictor values. If omitted,
#'      the original predictors are used.
#' @param type character indicating type of predictions: fitted means of response
#'      ("\code{response}"), corresponding linear predictor ("\code{link}"),
#'      fitted precision parameter phi ("\code{precision}"), fitted variances
#'      of response ("\code{variance}"), or fitted quantile(s) of the response
#'      distribution ("\code{quantile}").
#' @param at numeric vector indicating the level(s) at which quantiles should be
#'      predicted (only if \code{type = "quantile"}). Default is the median
#'      \code{at = 0.5}.
#' @param ... currently not used.
#'
#' @return Return a vector with the predicted values.
#'
#' @examples
#' \dontrun{
#' get(data("HIC", package = "robustbetareg"))
#' hic <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita | 1, data = HIC, alpha = 0.04)
#' cbind(predict(hic, type = "response"), predict(hic, type = "quantile", at = c(0.25, 0.5, 0.75)))}
#'
#' @export
predict = function(object, newdata = NULL, type = c("response", "link", "precision", "variance", "quantile"), at = 0.5, ...)
{
  UseMethod("predict")
}

#' @export
predict.robustbetareg = function(object, newdata = NULL, type = c("response", "link", "precision", "variance", "quantile"), at = 0.5, ...)
{
  type <- match.arg(type)
  if (type == "quantile") {
    qfun <- function(at, mu, phi) {
      rval <- sapply(at, function(p) stats::qbeta(p, mu * phi, (1 - mu) * phi))
      if (length(at) > 1L) {
        if (NCOL(rval) == 1L)
          rval <- matrix(rval, ncol = length(at), dimnames = list(unique(names(rval)),NULL))
        colnames(rval) <- paste("q_", at, sep = "")
      }
      else {
        rval <- drop(rval)
      }
      rval
    }
  }
  if (missing(newdata)) {
    rval <- switch(type, response = {
      object$fitted.values$mu.predict
    }, link = {
      set.link(object$link,object$link.phi)$linkfun.mu$linkfun(object$fitted.values$mu.predict)
    }, precision = {
      object$fitted.values$phi.predict
    }, variance = {
      object$fitted.values$mu.predict*(1-object$fitted.values$mu.predict)/(1+object$fitted.values$phi.predict)
    }, quantile = {
      qfun(at, object$fitted.values$mu.predict, object$fitted.values$phi.predict)
    })
    return(rval)
  }else{
    newdata=tryCatch(as.data.frame(newdata),error=function(e){newdata})
    x=model.matrix(object$formula,data=newdata,rhs = 1L)
    z=model.matrix(object$formula,data=newdata,rhs = 2L)

    rval <- switch (type, response = {
      set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean)
    }, link = {
      mu_predict=set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean)
      set.link(object$link,object$link.phi)$linkfun.mu$linkfun(mu_predict)
    }, precision = {
      set.link(object$link,object$link.phi)$linkfun.phi$inv.link(z%*%object$coefficients$precision)
    }, variance = {
      mu_predict=set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean)
      phi_predict=set.link(object$link,object$link.phi)$linkfun.phi$inv.link(z%*%object$coefficients$precision)
      mu_predict*(1-mu_predict)/phi_predict
    }, quantile={
      mu_predict=set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean)
      phi_predict=set.link(object$link,object$link.phi)$linkfun.phi$inv.link(z%*%object$coefficients$precision)
      qfun(at,mu_predict,phi_predict)
    })
    return(rval)
  }
}



#' @keywords internal
cooks.distance.robustbetareg=function(object,...)
{
  p=length(do.call("c",object$coefficients))
  linkobj=set.link(link.mu = object$link, link.phi = object$link.phi)
  y_hat=linkobj$linkfun.mu$inv.link(object$model$mean%*%object$coefficients$mean)
  MSE=as.numeric(t(object$y-y_hat)%*%(object$y-y_hat)/(object$n-p))
  D=NULL
  for(i in 1:object$n){
    fit.temp=robustbetareg(object$formula,data=object$data[-i,],alpha=object$Tuning)
    y_hat_temp=linkobj$linkfun.mu$inv.link(object$model$mean%*%fit.temp$coefficients$mean)
    D=c(D,t(y_hat-y_hat_temp)%*%(y_hat-y_hat_temp)/(MSE*p))
  }
  return(D)
}

#' @keywords internal
hatvalues.robustbetareg=function(object)
{
  mu_hat=object$fitted.values$mu.predict
  phi_hat=object$fitted.values$phi.predict
  y=object$y
  X=object$model$mean
  linkobj=set.link(link.mu=object$link,link.phi=object$link.phi)
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)

  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  W.PHI=diag(x=phi_hat*V_star*((d.link.mu)^(-2)))
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)
  return(diag(H))
}



