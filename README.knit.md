---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# robustbetareg

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/robustbetareg)](https://CRAN.R-project.org/package=robustbetareg)
[![Codecov test coverage](https://codecov.io/gh/yurimaluf/robustbetareg/branch/main/graph/badge.svg)](https://app.codecov.io/gh/yurimaluf/robustbetareg?branch=main)
<!-- badges: end -->

The **robustbetareg** package allows fitting robust beta regression. Currently,
four types of robust estimators are supported. They depend on a tuning constant 
which may be fixed or selected by a data-driven algorithm also implemented in the package.Diagnostic tools associated with the fitted model, such as the residuals and goodness-of-fit statistics, are implemented. Robust Wald-type tests are available.

## Installation

You can install the development version of **robustbetareg** from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yurimaluf/robustbetareg")
```

## Main function

The main function of the $\textbf{robustbetareg}$ package is <tt>robustbetareg()</tt>, which allows to fitting robust beta regression to proportional data; this explains the name. The arguments of <tt>robustbetareg()</tt> are:


```r
robustbetareg(formula, data, alpha, type = c("LSMLE", "LMDPDE", "SMLE", "MDPDE"),
  link = c("logit", "probit", "cloglog", "cauchit", "loglog"), link.phi = NULL,
  control = robustbetareg.control(...), model = TRUE, ... )
```

The <tt>robustbetareg()</tt> function returns an object of class "<tt>robustbetareg</tt>", similar to "<tt>betareg</tt>" and "<tt>glm</tt>" objects, for which some methods available. The <tt>summary()</tt> method returns a standard output, with coefficient estimates, standard errors, partial Wald tests and p values for the regression coefficients, the pseudo $R^2$, etc.. The <tt>type</tt> argument in <tt>robustbetareg()</tt> specifies the type of estimators to be used. The <tt>plot()</tt> method draws graphs for diagnostic analyses.

## Example


```r
library(robustbetareg)
## basic example code
```

In the following, an example is presented to illustrate the capacities of $\textbf{robustbetareg}$ package. We use the <tt>RiskManagerCost</tt> dataset, available in the package.
```
data("RiskManagerCost", package = "robustbetareg)
```

The response variable is <tt>FIRMCOST</tt> and the covariates are logarithm of total assets (<tt>SIZELOG</tt>) and a measure of the firm's industry risk (<tt>INDCOST</tt>). In the following, 
we fit the beta regression model using the maximum likelihood estimator and the LSMLE, a robust estimator with tuning constant selected by the data-driven algorithm.


```r
# MLE fit (fixed alpha equal to zero)
fit_MLE <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
                         data = RiskManagerCost, type = "LSMLE", alpha = 0,
                         link.phi = "log")
summary(fit_MLE)

# Choosing alpha via data-driven algorithm
fit_LSMLE <- robustbetareg(FIRMCOST ~ SIZELOG + INDCOST,
                            data = RiskManagerCost, type = "LSMLE",
                            link.phi = "log")
```

The goodness of fit is assessed using diagnostic graphs through the plot method.


```r
plot(fit_LSMLE)
```

Further details and examples on the R package $\textbf{robustbetareg}$ can be found using the help on R by typing:

```
help("robustbetareg")
```

## Reference

Maluf, Y.S., Ferrari, S.L.P., and Queiroz, F.F. (2022). Robust beta regression through the logit transformation. $\textit{arXiv}$:2209.11315.
