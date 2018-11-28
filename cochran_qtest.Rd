\name{cochran_qtest}
\alias{cochran_qtest}
\title{
cochran_qtest
}
\description{
The \code{cochran_qtest} function performs Cochran's Q-test for heterogeneity. Currently only two strata are supported, though more than two strata may be supported in the future.
}
\usage{
cochran_qtest(estimate1, lcl1, ucl1, estimate2, lcl2, ucl2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimate1}{
The ratio estimate (e.g. RR, OR, HR, IRR) from the first group/study.
}
  \item{lcl1}{
The lower confidence limit for the first group/study estimate.
}
  \item{ucl1}{
The upper confidence limit for the first group/study estimate.
}
  \item{measure2}{
The ratio estimate (e.g. RR, OR, HR, IRR) from the second group/study.
}
  \item{lcl2}{
The lower confidence limit for the second group/study estimate.
}
  \item{ucl2}{
The upper confidence limit for the second group/study estimate.
}
}
\details{
This function performs a statistical test for heterogeneity by calculating Cochran's Q-statistic. The null hypothesis is that any variation in estimates from the two studies is due to sampling variation, or minor differences in study methodology; in effect, the pooled estimate is equivalent to the stratum-specific estimates. The alternative hypothesis is that the study estimates are heterogeneous. This test is equivalent to the Wald test performed by the \code{wald.hgen.test} function also included in the \code{epimisc} package, and should produce an identical P-value when used to compare the same estimates.

The resulting statistic has a chi-square distribution with degrees of freedom equal to the number of strata - 1.
}
\value{
\item{pooled.estimate}{The pooled estimate}
\item{cochran.q}{Value of Cochran's Q-statistic}
\item{P}{P-value calculated using estimated Q-statistic}
}
\references{
Sedgwick P. Meta-analyses: tests of heterogeneity. BMJ. 2012; 344:e3971.
Kaufman JS, MacLehose RF. Which of these things is not like the others? Cancer. 2013; 119(24): 4216-4222.
}
\author{
Ryan Patrick Kyle
}
\examples{
## The function is currently defined as
function (estimate1, lcl1, ucl1, estimate2, lcl2, ucl2) 
{
    var1 <- (((log(ucl1) - log(estimate1))/qnorm(0.025)))^2
    var2 <- (((log(ucl2) - log(estimate2))/qnorm(0.025)))^2
    num_pool <- (log(estimate1)/var1) + (log(estimate2)/var2)
    denom_pool <- (1/var1) + (1/var2)
    pooled <- num_pool/denom_pool
    q <- (((log(estimate1) - pooled)^2)/var1) + (((log(estimate2) - pooled)^2)/var2)
    return(list(pooled.estimate = round(exp(pooled), 2), cochran.q = round(q, 
        2), P = round(1 - pchisq(q, 1), 4)))
  }
}
