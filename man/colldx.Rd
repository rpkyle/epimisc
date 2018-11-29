\name{colldx}
\alias{colldx}
\title{
Collinearity Diagnostics for Nonlinear Models Using the Information Matrix
}
\description{
The \code{colldx} function implements the Belsley, Kuh, and Welsch approach to diagnosing collinearity.
}
\usage{
colldx(mod, digits)
}
\arguments{
  \item{model}{
a model object for one of several supported model classes. 
}
  \item{digits}{
minimal number of significant digits, see \link{\code{print.default}}.
}
\details{
The \code{colldx} function is the R version of a SAS macro originally titled \code{collingenmodv9c}. Beau Bruce adapted the \code{colldiag} function from the \code{perturb} package to compute condition indices as in Matthew Zack's SAS macro.

The macro (and this R function) implement a method suggested by Belsley, Kuh, and Welsch in Regression Diagnostics: Identifying Influential Data and Sources of Collinearity (1980).

Several model classes are currently supported (though have not been extensively tested): \code{glm}, \code{coxph}, \code{geeglm}, \code{geese}, \code{gee}, \code{lmerMod}, \code{glmerMod}, \code{glmmML}, and \code{yagsResult}.
}
\author{
Written by Beau Bruce, modified for compatibility with various model classes by Ryan Patrick Kyle <ryan.kyle@mail.mcgill.ca>.
}
\references{
Belsley DA. A Guide to Using the Collinearity Diagnostics. Computer Science in Economics and Management 2014; 4:33-50.
}
\examples{
# As used by perturb/colldiag, example of SAS collinearity diagnostics from
# https://stats.idre.ucla.edu/sas/webbooks/reg/chapter2/regressionwith-saschapter-2-regression-diagnostics/
# 2.4 Tests for Collinearity
library(foreign)
elemapi <- read.dta("https://stats.idre.ucla.edu/stat/data/elemapi2.dta")
m2 <- lm(api00 ~ acs_k3+avg_ed+grad_sch+col_grad+some_col, data=elemapi)
colldx(m2, digits=3)
}