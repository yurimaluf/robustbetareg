% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustbetareg.R
\name{residuals.robustbetareg}
\alias{residuals.robustbetareg}
\alias{residuals}
\alias{.robustbetareg}
\title{Residuals Method for robustbetareg Objects}
\usage{
\method{residuals}{robustbetareg}(
  object,
  type = c("sweighted2", "pearson", "weighted", "sweighted", "sweighted.gamma",
    "sweighted2.gamma", "combined", "combined.projection"),
  ...
)
}
\arguments{
\item{object}{fitted model object of class \code{robustbetareg}.}

\item{type}{character indicating type of residuals to be used.}

\item{...}{currently not used.}
}
\value{
\code{residuals} returns a vector with the residuals of the type
     specified in the \code{type} argument.
}
\description{
The function provides several types of residuals for the robust beta regression
models: Pearson residuals (raw residuals scaled by square root of variance function)
and different kinds of weighted residuals proposed by Espinheira et al. (2008)
and Espinheira et al. (2017).
}
\details{
The definitions of the first four residuals are provided in
     Espinheira et al. (2008):  equation (2) for "\code{pearson}",
     equation (6) for "\code{weighted}", equation (7) for "\code{sweighted}",
     and equation (8) for "\code{sweighted2}". For the last four residuals
     the definitions are described in Espinheira et al. (2017): equations (7)
     and (10) for the "\code{sweighted.gamma}" and "\code{sweighted2.gamma}",
     respectively, equation (9) for "\code{combined}", and equation (11)
     for "\code{combined.projection}".
}
\examples{
\donttest{
get(data("HIC", package = "robustbetareg"))
fit.hic <- robustbetareg(HIC ~ URB + GDP | 1,
                         data = HIC, alpha = 0.04)
res <- residuals(fit.hic, type = "sweighted2")
#plot(res)
#abline(h = 0)
}

}
\references{
Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust
   beta regression through the logit transformation. \emph{Metrika}:61–81.\cr \cr
    Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). On Beta
    Regression Residuals. \emph{Journal of Applied Statistics}, 35:407–419.\cr \cr
    Espinheira, P.L., Santos, E.G.and Cribari-Neto, F. (2017). On nonlinear
    beta regression residuals. \emph{Biometrical Journal}, 59:445-461.\cr \cr
}
\seealso{
\code{\link[robustbetareg:robustbetareg]{robustbetareg}}
}
