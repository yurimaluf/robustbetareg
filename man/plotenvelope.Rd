% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustbetareg.R
\name{plotenvelope}
\alias{plotenvelope}
\title{Normal Probability Plots of Residuals with Simulated Envelope for robustbetareg Objects}
\usage{
plotenvelope(
  object,
  type = c("sweighted2", "pearson", "weighted", "sweighted", "sweighted.gamma",
    "sweighted2.gamma", "combined", "combined.projection"),
  conf = 0.95,
  n.sim = 100,
  PrgBar = TRUE,
  control = robustbetareg.control(...),
  ...
)
}
\arguments{
\item{object}{fitted model object of class \code{robustbetareg}.}

\item{type}{character indicating the type of residuals to be used, see
\code{\link{residuals.robustbetareg}}. Default is \code{type = "sweighted2"}.}

\item{conf}{numeric specifying the confidence level of the simulated
envelopes. Default is \code{conf = 0.95}.}

\item{n.sim}{a positive integer representing the number of iterations
to generate the simulated envelopes. Default is \code{n.sim = 100}.}

\item{PrgBar}{logical. If \code{PrgBar = TRUE} the progress bar will be shown
in the console. Default is \code{PrgBar = TRUE}.}

\item{control}{a list of control arguments specified via
\code{\link[robustbetareg:robustbetareg.control]{robustbetareg.control}}.}

\item{...}{arguments passed to \code{\link{plot}}.}
}
\value{
\code{plotenvelope} returns normal probability plot of residuals with simulated
     envelope.
}
\description{
\code{plotenvelope} is used to display normal probability plots of residuals
with simulated envelope for the robust beta regression. Currently, eight types of
residuals are supported: sweighted2, pearson, weighted, sweighted,
sweighted.gamma, sweighted2.gamma, combined, and combined.projection residuals.
}
\details{
The \code{plotenvelope} creates normal probability plots with simulated
     envelope (see Atkinson (1985) for details). Under the correct model,
     approximately 100*conf of the residuals are expected to be inside the
     envelope.
}
\examples{
\donttest{
get(data("HIC", package = "robustbetareg"))
hic <- robustbetareg(HIC ~ URB + GDP | GDP,
data = HIC, alpha = 0.06)
plotenvelope(hic, n.sim = 50)

get(data("Firm", package = "robustbetareg"))
rmc <- robustbetareg(FIRMCOST ~ INDCOST + SIZELOG | INDCOST + SIZELOG, data = Firm)
plotenvelope(rmc, conf = 0.90)}
}
\references{
Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust
   beta regression through the logit transformation. \emph{Metrika}:61–81.\cr \cr
    Atkinson, A.C. (1985) Plots, transformations and regression: an
    introduction to graphical methods of diagnostic regression analysis.
    \emph{Oxford Science Publications}, Oxford.
}
\seealso{
\code{\link{robustbetareg}}, \code{\link{robustbetareg.control}},
         \code{\link{residuals.robustbetareg}}
}
\author{
Yuri S. Maluf (\email{yurimaluf@gmail.com}),
Francisco F. Queiroz (\email{ffelipeq@outlook.com}) and Silvia L. P. Ferrari.
}
