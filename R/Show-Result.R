#' @export
print.robustbetareg=function(x,...)
{
  object=x
  cat("Call: \n")
  print(object$call)
  cat("\n")
  cat("Coefficients (mean model with",object$link,"link):\n")
  print(object$coefficients$mean)
  cat("\n")
  cat("Coefficients (precision model with",object$link.phi,"link):\n")
  print(object$coefficients$precision)
  if(!object$converged)
  {
    cat("\n")
    cat("The algorithm did not reach the convergence.\n")
  }
  cat("------\n")
  cat(paste0("Tuning value: alpha=",object$Tuning))
}


#' @export
print.WaldTest_robustbetareg=function(x, digits = max(3, getOption("digits") - 3),...)
{
  object=x
  if(object$general)
  {
    cat("-- Wald Type Test -- \n")
    if(!is.null(object$beta.wald))
    {
      p.valor=round(object$beta.wald$pValue, digits = digits)
      obs=star.obs(p.valor)
      if(p.valor<=2e-16){p.valor="<2e-16"}
      if(p.valor>2e-16){p.valor=paste0("=",p.valor)}

      cat("Null Hypothesis: all mean coefficients equal to zero \n")
      cat(paste0("Value=",formatC(object$beta.wald$W.alpha, digits = digits),", df=",object$beta.wald$df,", p-Value",p.valor,obs,"\n"))
    }
    if(!is.null(object$gamma.wald))
    {
      p.valor=round(object$gamma.wald$pValue, digits = digits)
      obs=star.obs(p.valor)
      if(p.valor<=2e-16){p.valor="<2e-16"}
      if(p.valor>2e-16){p.valor=paste0("=",object$gamma.wald$pValue)}

      cat("Null Hypothesis: all precision coefficients equal to zero \n")
      cat(paste0("Wald test=",formatC(object$gamma.wald$W.alpha, digits = digits),", df=",object$gamma.wald$df,", p-Value",p.valor,obs,"\n"))
    }
  }else{
    cat("-- Wald Type Test -- \n")
    p.valor=round(object$pValue, digits = digits)
    obs=star.obs(p.valor)
    if(p.valor>0.0001){p.valor=paste0("=",p.valor)}
    else if(p.valor<=2e-16){p.valor="<2e-16"}

    cat("Null Hypothesis: set by the user \n")
    cat(paste0("Wald test=",formatC(object$W.alpha, digits = digits),", df=",object$df,", p-Value",p.valor,obs,"\n"))
  }
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  cat(paste0("Results based on ",object$method, "\n"))
}


#' Methods for robustbetareg Objects
#'
#' Methods for extracting information from fitted robust beta regression model objects of class "\code{robustbetareg}"
#'
#' @param object fitted model of class \code{robustbetareg}.
#' @param type character specifying type of residuals to be included in the summary output, see \code{\link[robustbetareg:residuals]{residuals.robustbetareg}}.
#' @param ... currently not used.
#'
#' @details A set of standard extractor functions for fitted model objects is available for objects of class "\code{robustbetareg}",
#' including methods to the generic functions print and summary which print the estimated coefficients along with some further information.
#'
#' @references \href{https://www.tandfonline.com/doi/abs/10.1080/02664760701834931}{Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). On Beta Regression Residuals. Journal of Applied Statistics, 35(4), 407–419.}
#' @references \href{https://www.tandfonline.com/doi/abs/10.1080/02664760701834931}{Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). On Beta Regression Residuals. Journal of Applied Statistics, 35(4), 407–419.}
#' @references \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201600136}{Espinheira, P.L., Santos, E.G.and Cribari-Neto, F. (2017). On nonlinear beta regression residuals. Biometrical Journal, 59(3), 445-461.}
#' @references \href{https://doi.org/10.1177/0962280217738142}{Ghosh, A. Robust inference under the beta regression model with application to health care studies. Statistical Methods in Medical Research, 28:271-888 (2019).}
#' @references \href{https://doi.org/10.48550/arXiv.2209.11315}{Maluf, Y. S., Ferrari, S. L., & Queiroz, F. F. (2022). Robust beta regression through the logit transformation. arXiv}
#' @references \href{https://doi.org/10.1007/s00362-022-01320-0}{Ribeiro, K. A. T. Ferrari, S. L. P. Robust estimation in beta regression via maximum Lq-likelihood. Statistical Papers (2022).}
#'
#' @seealso \code{\link[robustbetareg:robustbetareg]{robustbetareg}}
#'
#' @examples
#' \dontrun{data("HIC", package="robustbetareg")
#' fit=robustbetareg(Percent_HIC~Urbanization+GDP_percapita|1,data=HIC,alpha=0.06)
#' summary(fit)
#' coef(fit)
#' }
#'
#'@export
summary.robustbetareg=function(object, type = "sweighted2", ...)
{
  type <- match.arg(type, c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3"))
  object$residuals = residuals(object,type=type)
  object$residuals.type <- type
  k <- length(object$coefficients$mean)
  m <- length(object$coefficients$precision)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], precision = cf[seq.int(length.out = m) + k, , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$precision) <- names(object$coefficients$precision)
  object$coefficients <- cf
  class(object) <- "summary.robustbetareg"
  object
}

#' @rdname summary.robustbetareg
#'
#' @param model character specifying for which component of the model coefficients/covariance should be extracted.
#'
#' @export
coef.robustbetareg=function(object,model=c("full","mean","precision"),...)
{
  cf <- object$coefficients
  model=match.arg(model)
  switch(model, mean = {
    cf$mean
  }, precision = {
    cf$precision
  }, full = {
    nam1 <- names(cf$mean)
    nam2 <- names(cf$precision)
    cf <- c(cf$mean, cf$precision)
    names(cf) <- c(nam1, if (identical(nam2, "(Phi)")) "(phi)" else paste("(phi)",nam2, sep = "_"))
    cf
  })
}

#' @export
#' @importFrom stats printCoefmat quantile
print.summary.robustbetareg=function(x, digits = max(3, getOption("digits") - 3), ...){
  #browser()
  cat("Call:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if (!x$converged) {
    cat("model did not converge\n")
  }else {
    types <- c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3")
    Types <- c("Standard Weighted 2","Pearson","Weighted","Sweighted","Sweighted Gamma","Sweighted2 Gamma","Combined","Combined Projection","Standard Weighted 3")
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(stats::quantile(x$residuals)), digits = digits), .Names = c("Min", "1Q", "Median","3Q", "Max")))
    if (NROW(x$coefficients$mean)) {
      cat(paste("\nCoefficients (mean model with ", x$link, " link):\n", sep = ""))
      stats::printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
    }
    else cat("\nNo coefficients (in mean model)\n")
    if (NROW(x$coefficients$precision)) {
        cat(paste("\nPhi coefficients (precision model with ", x$link.phi, " link):\n", sep = ""))
        printCoefmat(x$coefficients$precision, digits = digits, signif.legend = FALSE)
    }else cat("\nNo coefficients (in precision model)\n")
    if (getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    if(x$Tuning==0){
      cat("\nType of estimator:", "MLE")
    }else{
      cat("\nType of estimator:", switch(x$method, LSMLE = "LSMLE",  LMDPDE = "LMDPDE", SMLE = "SMLE",  MDPDE = "MDPDE"))
    }
    cat(paste0("\nTuning value: alpha=", x$Tuning))
    if(x$Optimal.Tuning){
      cat("\nTuning of", x$method ,"generated by the data-driven algorithm.")
      if(!is.null(x$message))
      {
        if(x$message=="Lack of stability"){cat("\nLack of stability.")}
      }
    }
    if (!is.na(x$pseudo.r.squared))
      cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
  }
  invisible(x)
}

#' Interactive plots for diagnostic of robust betareg models
#'
#' Several types of standard diagnostic plots can be produced interactively, involving various kinds of residuals, influence measures, weights etc.
#'
#' @param x fitted model object of class \code{robustbetareg}.
#' @param ask logical. If TRUE the user is asked before each plot.
#' @param ... other parameters to be passed through to plotting functions.
#'
#' @examples
#' \dontrun{
#' get(data("HIC", package = "robustbetareg"))
#' hic <- robustbetareg(Percent_HIC ~ Urbanization + GDP_percapita | GDP_percapita,
#' data = HIC, alpha = 0.06)
#' plot(hic)}
#'
#' @importFrom stats residuals
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics identify
#'
#' @export
plot.robustbetareg=function(x,ask=TRUE,...)
{
  object=x
  getinfo=Sys.info()
  user=getinfo[which(names(getinfo)=="user")]
  text.main2="the graph number >\n [1] Residuals \n [2] Residuals x Linear predictor \n [3] Cook's Distance \n [4] Weights \n [5] Weigths x Residuals \n [0] Exit \n"
  text.main=paste("Select",text.main2)
  if(!is.na(user))
  {
    user=paste0(toupper(substring(user,1,1)),substring(user,2))
    text.main=paste0(user,", select ",text.main2)
  }
  text.n1="Select the residual number: \n [1] sweighted2 \n [2] sweighted \n [3] pearson \n [4] weighted \n [5] sweighted.gamma \n [6] sweighted2.gamma \n [7] combined \n [8] combined.projection \n [9] back \n [0] exit \n"
  text.n2="Select the residual type to match with linear predictor: \n [1] sweighted2 \n [2] sweighted \n [3] pearson \n [4] weighted \n [5] sweighted.gamma \n [6] sweighted2.gamma \n [7] combined \n [8] combined.projection \n [9] back \n [0] exit \n"
  text.n5="Select the residual type to match with weights: \n [1] sweighted2 \n [2] sweighted \n [3] pearson \n [4] weighted \n [5] sweighted.gamma \n [6] sweighted2.gamma \n [7] combined \n [8] combined.projection \n [9] back \n [0] exit \n"
  show=TRUE
  while(show){
    show.1=show.2=show.5=TRUE
    rstudioapi::sendToConsole("",execute=F,focus=T,echo=T)
    n <- as.numeric(readline(cat(crayon::green(text.main))))
    if(n==1)
    {
      while(show & show.1)
      {
        m <-readline(prompt = cat(crayon::green(text.n1)))
        if(m==1)
        {
          res=residuals(object,type="sweighted2")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted 2 Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==2)
        {
          res=residuals(object,type="sweighted")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==3)
        {
          res=residuals(object,type="pearson")
          plot(res,xlab="Obs. number",ylab="Pearson Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==4)
        {
          res=residuals(object,type="weighted")
          plot(res,xlab="Obs. number",ylab="Weighted Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==5)
        {
          res=residuals(object,type="sweighted.gamma")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted Gamma Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==6)
        {
          residuals(object,type="sweighted2.gamma")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted 2 Gamma Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==7)
        {
          res=residuals(object,type="combined")
          plot(res,xlab="Obs. number",ylab="Combined Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==8)
        {
          res=residuals(object,type="combined.projection")
          plot(res,xlab="Obs. number",ylab="Combined Projection Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==0){show=FALSE}
        if(m==9){show.1=FALSE}
      }
    }
    if(n==2)
    {
      while(show & show.2)
      {
        m <-readline(prompt = cat(crayon::green(text.n2)))
        if(m==1)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted2"),xlab="Linear predictor",ylab="Standardized Weighted 2 Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==2)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted"),xlab="Linear predictor",ylab="Standardized Weighted Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==3)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="pearson"),xlab="Linear predictor",ylab="Pearson Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==4)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="weighted"),xlab="Linear predictor",ylab="Weighted Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==5)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted.gamma"),xlab="Linear predictor",ylab="Standardized Weighted Gamma Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==6)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted2.gamma"),xlab="Linear predictor",ylab="Standardized Weighted 2 Gamma Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==7)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="combined"),xlab="Linear predictor",ylab="Combined Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==8)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="combined.projection"),xlab="Linear predictor",ylab="Combined Projection Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==0){show=FALSE}
        if(m==9){show.2=FALSE}
      }
    }
    if(n==3)
    {
      plot(stats::cooks.distance(object),type="h",xlab="Obs. number",ylab="Cook's distance",main="Cook's distance plot",...)
    }
    if(n==4)
    {
      plot(object$weights,xlab="Obs. number",ylab="Weights",main = "Weights plot",...)
      abline(h=0)
    }
    if(n==5)
    {
      while(show & show.5)
      {
        m <-readline(prompt = cat(crayon::green(text.n5)))
        if(m==1)
        {
          plot(residuals(object,type="sweighted2"),object$weights,ylab="Weights",xlab="Standardized Weighted 2 Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==2)
        {
          plot(residuals(object,type="sweighted"),object$weights,ylab="Weights",xlab="Standardized Weighted Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==3)
        {
          plot(residuals(object,type="pearson"),object$weights,ylab="Weights",xlab="Pearson Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==4)
        {
          plot(residuals(object,type="weighted"),object$weights,ylab="Weights",xlab="Weighted Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==5)
        {
          plot(residuals(object,type="sweighted.gamma"),object$weights,ylab="Weights",xlab="Standardized Weighted Gamma Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==6)
        {
          plot(residuals(object,type="sweighted2.gamma"),object$weights,ylab="Weights",xlab="Standardized Weighted 2 Gamma Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==7)
        {
          plot(residuals(object,type="combined"),object$weights,ylab="Weights",xlab="Combined Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==8)
        {
          plot(residuals(object,type="combined.projection"),object$weights,ylab="Weights",xlab="Combined Projection Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==0){show=FALSE}
        if(m==9){show.5=FALSE}
      }
    }
    if(n==0)
    {
      show=FALSE
    }
  }
}


