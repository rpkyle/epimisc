\name{std_measure}
\alias{std_measure}
\title{
Standardization of epidemiologic measures of association
}
\description{
This function may be used to standardize risk ratios and risk differences by one of several methods. Odds ratios and rate ratios are not yet supported, but will likely be valid options for the \code{measure} statement in a future release. This function was inspired by code provided online by Bernhard Klingenberg, and cited below in the references.
}
\usage{
std_measure(dataset, conflev = 0.95, method, measure, user.weights = NA, digits)
}
\arguments{
  \item{dataset}{
A dataset in tabular form, with one or more strata, as produced by the \code{table} function. The rows should represent exposure status, starting with the referent group; the columns should represent levels of the outcome variable, also starting with the referent.
}
  \item{conflev}{
The desired confidence level. The default is to estimate 95\% confidence intervals.
}
  \item{method}{
The estimation method for the standardized result; valid choices are "MH" for Mantel-Haenszel standardization, "internal" for standardization to the exposed, "external" for standardization to the unexposed, or 'user.weights' to use a vector of stratum-specific weights provided by the user.
}
  \item{measure}{
The desired measure to estimate; valid options are currently limited to "risk.ratio" and "risk.difference".
}
  \item{user.weights}{
A vector of user-supplied weights with length equal to the number of strata.
}
  \item{digits}{
The number of decimal places to report in output. The default is 3.
}
}
\references{
Klingenberg B. A new and improved confidence interval for the Mantel-Haenszel risk difference. Stat Med. 2014; 33(17):2968-2983.

Rothman KJ, Greenland S, Lash TL. Modern Epidemiology. 3rd ed. Philadelphia, PA: Lippincott, Williams & Wilkins; 2012.
}
\author{
Ryan Patrick Kyle <ryan.kyle@mail.mcgill.ca>, Bernhard Klingenberg (authored \code{stratMHRD.R} function, upon which \code{std_measure} is heavily based).
}
\examples{
# borrowing a Stata dataset to use for example
require(foreign)
ugdp <- read.dta("http://www.stata-press.com/data/r11/ugdp.dta")

# crude trick to "expand" a frequency-weighted table for use
# with epi.2by2 and this function
rows <- rep(1:nrow(ugdp), ugdp[, 4])
ugdp <- ugdp[rows,][,-4]

# get the data into proper shape for use with function
# note that this is the same configuration required by epi.2by2
# from epiR package:
ugdp$case <- factor(ugdp$case, levels = c(1, 0), labels = c("case", "non-case"))
names(ugdp)[2] <- "outcome status"
ugdp$exposed <- factor(ugdp$exposed, levels = c(1, 0), labels = c("exposed", "unexposed"))
names(ugdp)[3] <- "exposure status"

tab1 <- table(ugdp[,"exposure status"], ugdp[,"outcome status"], ugdp[,"age"], 
   dnn = c("Exposure status", "Outcome status", "Age"))
print(tab1)

# weighted risk ratio, to distribution among the exposed (internal standardization)
# cf. "cs case exposed [fw=pop], by(age) istandard" in Stata
std_measure(tab1, measure = "risk.ratio", method = "internal", digits = 6)

# weighted risk difference, to distribution among the unexposed (external standardization)
# cf. "cs case exposed [fw=pop], by(age) rd estandard" in Stata
std_measure(tab1, measure = "risk.difference", method = "external", digits = 6)

# user-specified weights, equivalent to above, weighted risk difference
std_measure(tab1, measure = "risk.difference", method = "user.weights", user.weights = c(120, 85), digits = 6)

# Risk ratio using Mantel-Haenszel approach
# cf. "cs case exposed [fw=pop], by(age)" 
std_measure(tab1, measure = "risk.ratio", method = "MH", digits = 6)
}
