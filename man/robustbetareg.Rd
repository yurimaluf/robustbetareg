% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustbetareg.R
\name{robustbetareg}
\alias{robustbetareg}
\alias{LMDPDE.fit}
\alias{LSMLE.fit}
\alias{MDPDE.fit}
\alias{SMLE.fit}
\title{Robust Beta Regression}
\usage{
robustbetareg(
  formula,
  data,
  alpha,
  type = c("LSMLE", "LMDPDE", "SMLE", "MDPDE"),
  link = c("logit", "probit", "cloglog", "cauchit", "loglog"),
  link.phi = NULL,
  control = robustbetareg.control(...),
  model = TRUE,
  ...
)

LMDPDE.fit(y, x, z, alpha = NULL, link = "logit",
link.phi = "log", control = robustbetareg.control(...), ...)

LSMLE.fit(y, x, z, alpha = NULL, link = "logit",
link.phi = "log", control = robustbetareg.control(...), ...)

MDPDE.fit(y, x, z, alpha = NULL, link = "logit",
link.phi = "log", control = robustbetareg.control(...), ...)

SMLE.fit(y, x, z, alpha = NULL, link = "logit",
link.phi = "log", control = robustbetareg.control(...), ...)
}
\arguments{
\item{formula}{symbolic description of the model. See Details for further
information.}

\item{data}{dataset to be used.}

\item{alpha}{numeric in \eqn{[0,1)} indicating the value of the tuning constant
alpha. \code{alpha = 0} leads to the maximum likelihood estimator.
Robust procedures require \code{alpha} greater than zero.
If this argument is suppressed, the tuning constant will be selected
automatically through the data-driven algorithm proposed by Ribeiro and
Ferrari (2022).}

\item{type}{character specifying the type of robust estimator to be used in the
estimation process. Supported estimators are "\code{LSMLE}" ,
"\code{LMDPDE}",  "\code{SMLE}", and "\code{MDPDE}"; for details, see Maluf
 et al. (2022). The "\code{LSMLE}" is the default.}

\item{link}{an optional character that specifies the link function of the
mean submodel (mu). The "\code{logit}", "\code{probit}", "\code{cloglog}",
"\code{cauchit}", "\code{loglog}" functions are supported. The \code{logit}
function is the default.}

\item{link.phi}{an optional character that specifies the link function of the
precision submodel (phi). The "\code{identity}", "\code{log}", "\code{sqrt}"
functions are supported. The default is \code{log} unless formula is of type
\code{y ~ x} where the default is "\code{identity}".}

\item{control}{a list of control arguments specified via
\code{\link{robustbetareg.control}}.}

\item{model}{logical. If \code{TRUE} the corresponding components of the fit
(model frame, response, model matrix) are returned.}

\item{...}{argument to be passed to \code{\link{robustbetareg.control}}.}

\item{y, x, z}{\code{y} must be a numeric response vector (with values in
\eqn{(0,1)}), \code{x} must be a numeric regressor matrix for the mean
 submodel, and \code{z} must be a numeric regressor matrix for the precision
 submodel.}
}
\value{
\code{robustbetareg} returns an object of class "\code{robustbetareg}" with a list of the following components:\tabular{ll}{
   \code{coefficients} \tab a list with the "\code{mean}" and "\code{precision}"
   coefficients. \cr
   \tab \cr
   \code{vcov} \tab covariance matrix. \cr
   \tab \cr
   \code{converged} \tab  logical indicating successful convergence of the
      iterative process. \cr
   \tab \cr
   \code{fitted.values} \tab a vector with the fitted values of the mean submodel. \cr
   \tab \cr
   \code{start} \tab a vector with the starting values used in the iterative process. \cr
   \tab \cr
   \code{weights} \tab the weights of each observation in the estimation process. \cr
   \tab \cr
   \code{Tuning} \tab value of the tuning constant (automatically chosen or fixed) used
      in the estimation process. \cr
   \tab \cr
   \code{residuals} \tab a vector of standardized weighted residual 2 (see Espinheira et al. (2008)). \cr
   \tab \cr
   \code{n} \tab number of observations. \cr
   \tab \cr
   \code{link} \tab link function used in the mean submodel. \cr
   \tab \cr
   \code{link.phi} \tab link function used in the precision submodel. \cr
   \tab \cr
   \code{Optimal.Tuning} \tab logical indicating whether the data-driven algorithm
      was used. \cr
   \tab \cr
   \code{pseudo.r.squared} \tab pseudo R-squared value. \cr
   \tab \cr
   \code{control} \tab the control arguments passed to the data-driven algorithm and
     \code{optim} call. \cr
   \tab \cr
   \code{std.error} \tab the standard errors. \cr
   \tab \cr
   \code{method} \tab type of estimator used. \cr
   \tab \cr
   \code{call} \tab the original function call. \cr
   \tab \cr
   \code{formula} \tab the formula used. \cr
   \tab \cr
   \code{model} \tab the full model frame. \cr
   \tab \cr
   \code{terms} \tab a list with elements "\code{mean}", "\code{precision}" and "\code{full}"
       containing the term objects for the respective models.  \cr
   \tab \cr
   \code{y} \tab the response variable. \cr
   \tab \cr
   \code{data} \tab the dataset used. \cr
}
}
\description{
Fit robust beta regression models for rates and proportions via LSMLE, LMDPDE,
 SMLE and MDPDE. Both mean and precision of the response variable are modeled
 through parametric functions.
}
\details{
Beta regression models are employed to model continuous response
   variables in the unit interval, like rates and proportions. The maximum
   likelihood-based inference suffers from
   the lack of robustness in the presence of outliers. Based on
   the density power divergence, Ghosh (2019) proposed the minimum density
   power divergence estimator (MDPDE). Ribeiro and Ferrari (2022) proposed an
   estimator based on the maximization of a reparameterized Lq-likelihood;
   it is called SMLE. These estimators require suitable restrictions in the
   parameter space. Maluf et al. (2022) proposed robust estimators based on
   the MDPDE and the SMLE which have the advantage of overcoming this drawback.
   These estimators are called LMDPDE and LSMLE. For details, see the
   cited works. The four estimators are implemented in the \code{robustbetareg}
   function. They depend on a tuning constant (called \eqn{\alpha}).
   When the tuning constant is fixed and equal to 0, all of the estimators
   coincide with the maximum likelihood estimator. Ribeiro and Ferrari (2022)
   and Maluf et al. (2022) suggest using a data-driven algorithm to select the
   optimum value of \eqn{\alpha}. This algorithm is implemented in
    \code{robustbetareg} by default when the argument "\code{alpha}" is
    suppressed.\cr \cr
    The formulation of the model has the same structure as in the usual functions
    \code{\link[stats]{glm}} and \code{\link[betareg]{betareg}}. The argument
    \code{formula} can comprise of three parts (separated by the symbols
    "\eqn{~}" and "\eqn{|}"), namely: observed response variable in the unit
    interval, predictor of the mean submodel, with link function \code{link}
    and predictor of the precision submodel, with \code{link.phi}
    link function. If the model has constant precision, the third part may be
    omitted and the link function for phi is "\code{identity}" by default.
    The tuning constant \code{alpha} may be treated as fixed or not (chosen
    by the data-driven algorithm). If \code{alpha} is fixed, its value
    must be specified in the \code{alpha} argument. \cr \cr
    Some methods are available for objects of class "\code{robustbetareg}",
    see \code{\link{plot.robustbetareg}}, \code{\link{summary.robustbetareg}},
    \code{\link{coef.robustbetareg}}, and \code{\link{residuals.robustbetareg}},
     for details and other methods.
}
\examples{
#### Risk Manager Cost data
data("Firm")

# MLE fit (fixed alpha equal to zero)
fit_MLE <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
                         data = Firm, type = "LMDPDE", alpha = 0)
summary(fit_MLE)
\donttest{
# MDPDE with alpha = 0.04
fit_MDPDE <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
                           data = Firm, type = "MDPDE",
                           alpha = 0.04)
summary(fit_MDPDE)

# Choosing alpha via data-driven algorithm
fit_MDPDE2 <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
                            data = Firm, type = "MDPDE")
summary(fit_MDPDE2)

# Similar result for the LMDPDE fit:
fit_LMDPDE2 <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
                             data = Firm, type = "LMDPDE")
summary(fit_LMDPDE2)

# Diagnostic plots


#### HIC data
data("HIC")

# MLE (fixed alpha equal to zero)
fit_MLE <- robustbetareg(HIC ~ URB + GDP |
                         GDP, data = HIC, type = "LMDPDE",
                         alpha = 0)
summary(fit_MLE)

# SMLE and MDPDE with alpha selected via data-driven algorithm
fit_SMLE <- robustbetareg(HIC ~ URB + GDP |
                          GDP, data = HIC, type = "SMLE")
summary(fit_SMLE)
fit_MDPDE <- robustbetareg(HIC ~ URB + GDP |
                           GDP, data = HIC, type = "MDPDE")
summary(fit_MDPDE)
# SMLE and MDPDE return MLE because of the lack of stability

# LSMLE and LMDPDE with alpha selected via data-driven algorithm
fit_LSMLE <- robustbetareg(HIC ~ URB + GDP |
                           GDP, data = HIC, type = "LSMLE")
summary(fit_LSMLE)
fit_LMDPDE <- robustbetareg(HIC ~ URB + GDP |
                            GDP, data = HIC, type = "LMDPDE")
summary(fit_LMDPDE)
# LSMLE and LMDPDE return robust estimates with alpha = 0.06


# Plotting the weights against the residuals - LSMLE fit.
plot(fit_LSMLE$residuals, fit_LSMLE$weights, pch = "+", xlab = "Residuals",
    ylab = "Weights")

# Excluding outlier observation.
fit_LSMLEwo1 <- robustbetareg(HIC ~ URB + GDP |
                              GDP, data = HIC[-1,], type = "LSMLE")
summary(fit_LSMLEwo1)

# Normal probability plot with simulated envelope
 plotenvelope(fit_LSMLE)
}
}
\references{
Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust
   beta regression through the logit transformation. \emph{Metrika}:61–81.\cr \cr
   Ribeiro, T.K.A. and Ferrari, S.L.P.  (2022). Robust estimation in beta regression
   via maximum Lq-likelihood. \emph{Statistical Papers}. DOI: 10.1007/s00362-022-01320-0. \cr \cr
   Ghosh, A. (2019). Robust inference under the beta regression model with
   application to health care studies. \emph{Statistical Methods in Medical
   Research}, 28:271-888.\cr \cr
   Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). On beta regression residuals. \emph{Journal of Applied Statistics}, 35:407–419.
}
\seealso{
\code{\link{robustbetareg.control}}, \code{\link{summary.robustbetareg}}, \code{\link{residuals.robustbetareg}}
}
\author{
Yuri S. Maluf (\email{yurimaluf@gmail.com}), Francisco F. Queiroz (\email{ffelipeq@outlook.com}) and Silvia L. P. Ferrari.
}
