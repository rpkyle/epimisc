\name{lincomR}
\alias{lincomR}
\title{
Linear combinations of estimators
}
\description{
The \code{lincomR} function is similar to the Stata \code{lincom} function. This function calculates point estimates and their standard errors for linear combinations of coefficients after fitting a model. Exponentiated estimates and their confidence intervals are also provided when appropriate.
}
\usage{
lincomR(model, contrast, conflev = 0.95, digits = 3)
}
\arguments{
  \item{model}{
A model object. Currently the following classes are supported: \code{coxph}, \code{glm}, \code{geese}, \code{geeglm}, \code{glmmML}, \code{mer}, and \code{yagsResult}.
}
  \item{contrast}{
A contrast vector of integer values, equal to the number of parameters estimated by the model, including the intercept (except for Cox proportional hazards models).
}
  \item{conflev}{
The confidence level used for estimating confidence intervals. The default is to estimate 95\% confidence intervals.
}
  \item{digits}{
The number of decimal places to report for estimates. The default is 3.
}
}
\details{
This function takes four inputs: a model from a Cox proportional hazards, generalized linear, generalized linear mixed, or generalized linear model fit with GEE, a contrast vector, a confidence level, e.g. 0.95 for 95\% confidence intervals, and (optionally) the number of decimal places to report when displaying estimates. The contrast vector should be of length = number of model covariates + 1 (for the intercept -- except for Cox proportional hazards models), and this contrast represents the difference between the two specified models. Values in the contrast vector should be negative if the parameter is to be subtracted, and the value should represent any constant to be multiplied by (for continuous variables). Function returns estimate, its standard error, Wald z-score and associated P-value, and confidence interval.                                                  
}
\value{
\item{conflev}{The specified confidence level for the confidence intervals estimated by the function. The default is to estimate 95\% confidence intervals.}
\item{est}{The point estimate for the chosen contrast.}
\item{se.est}{The standard error of the point estimate prior to exponentiation.}
\item{wald.z}{Z-score corresponding to Wald test.}
\item{pvalue}{P-value corresponding to Wald test.}
\item{exp.est}{The exponentiated (ratio) point estimate, if link is \code{log} or \code{logit}, or model type is Cox PH.}
\item{se.exp.est}{The standard error of the exponentiated estimate.}
\item{lower.bound}{The lower bound of the confidence interval for \code{exp.est} corresponding to \code{conflev}. If link function is not \code{log}, \code{logit}, or model type was not \code{coxph}, the bounds will not be exponentiated.}
\item{upper.bound}{The upper bound of the confidence interval for \code{exp.est} corresponding to \code{conflev}. See note above for \code{lower.bound}.}
}
\author{
Originally written by Kathryn Morrison <kt.morrison@mail.mcgill.ca>, modified for compatibility with various model classes and to provide nicely formatted output by Ryan Patrick Kyle <ryan.kyle@mail.mcgill.ca>.
}
\note{
This function is currently *experimental* and is not guaranteed to work properly with all model classes and link functions, and has not been extensively tested with mixed model objects. User feedback is appreciated.
}

\examples{
# borrowing the same example used in Stata 13 documentation for --lincom--
# http://www.stata.com/manuals13/rlincom.pdf
require(foreign)
lbw3 <- read.dta("http://www.stata-press.com/data/r11/lbw3.dta")

# logit low age lwd i.race smoke ptd ht ui
lbwmod <- glm(low ~ age + lwd + as.factor(race) + smoke + ptd + ht + ui, 
	family = binomial(link="logit"), data = lbw3)

# odds ratio for smoking among race == 1 (black), relative to non-smoking among race == 0 (white)
lincomR(lbwmod, c(0, 0, 0, 1, 0, 1, 0, 0, 0), conflev = 0.95, digits = 6)
}
