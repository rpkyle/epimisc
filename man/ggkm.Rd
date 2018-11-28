\name{ggkm}
\alias{ggkm}
\title{
Generate Kaplan-Meier plot using ggplot2 and gridExtra
}
\description{
Given a \code{survfit} object produced by the \code{survfit} function, \code{ggkm} will produce an attractive Kaplan-Meier plot, with an optional table summarizing the number of individuals at risk along the footer of the plot.
}
\usage{
ggkm(sfit, table = FALSE, xlabs = "Time", ylabs = "survival probability", 
ystratalabs = NULL, ystrataname = NULL, timeby = NULL, main = "Kaplan-Meier Plot", 
pval = TRUE)
}
\arguments{
  \item{sfit}{
  An object produced by the \code{survfit} function.
}
  \item{table}{
Boolean value; if true, a table with number at risk is appended to the bottom of the plot. Default is \code{FALSE}.
}
  \item{xlabs}{
 Character string used to label x-axis.
  }
  \item{ylabs}{
 Character string used to label y-axis.
  }
  \item{ystrataname}{
 Character string used to name strata legend; if not supplied, "Strata" is used instead.
  }
  \item{timeby}{
 Specifies interval used to define breakpoints for at-risk table.
  }
  \item{main}{
 Character string used to title finished plot; if not specified, default is "Kaplan-Meier Plot".
  }
  \item{pval}{
 Boolean value; if true, the P-value for the chi-square (log-rank) test of equality performed by \code{survdiff} is printed on the finished plot. Default is \code{TRUE}.
  }
}
\value{
Produces a plot created with \code{ggplot2}.
}
\author{
Abhijit Dasgupta,
Gary Collins,
Gil Tomas,
Mark Cowley,
Dieter Menne
}
\source{
\code{https://statbandit.wordpress.com/2011/03/08/an-enhanced-kaplan-meier-plot/}
 
\code{https://statbandit.wordpress.com/2014/04/01/kaplan-meier-plots-using-ggplots2-updated/}

\code{https://github.com/dmenne/dmisc2/tree/master}
}

\examples{
require(survival)
mfit <- survfit(Surv(futime, death) ~ sex, data=mgus)

# no table of individuals at-risk
ggkm(mfit)

# with table of individuals at-risk
ggkm(mfitm, table = TRUE)
}
