ggkm <- function(sfit, table = FALSE,
                      xlabs = "Time", ylabs = "survival probability",
                      ystratalabs = NULL, ystrataname = NULL,
                      timeby = NULL, main = "Kaplan-Meier Plot",
                      pval = TRUE) {
  surv <- NULL
  n.risk <- NULL
  if (is.null(ystratalabs)) {
    ystratalabs <- as.character(levels(summary(sfit)$strata))
  }
  m <- max(nchar(ystratalabs))
  if (is.null(ystrataname))
    ystrataname <- "Strata"
  if (is.null(timeby)) {
    times <- pretty(sfit$time)
  } else {
    times <- seq(0, max(sfit$time), by = timeby)
  }
  .df <- data.frame(
    time = sfit$time, n.risk = sfit$n.risk,
    n.event = sfit$n.event, surv = sfit$surv, strata = summary(sfit, censored = T)$strata,
    upper = sfit$upper, lower = sfit$lower
  )
  levels(.df$strata) <- ystratalabs
  zeros <-
    data.frame(
      time = 0, surv = 1, strata = factor(ystratalabs, levels = levels(.df$strata)),
      upper = 1, lower = 1
    )
  .df <- rbind.fill(zeros, .df) # This should be replaced, requires import of plyr
  d <- length(levels(.df$strata))
  p <- ggplot(.df, aes(time, surv, group = strata)) +
    geom_step(aes(linetype = strata), size = 0.7) +
    theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    xlab(xlabs) +  ylab(ylabs) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position = c(ifelse(m < 10, .28, .35), ifelse(d < 4, .25, .35))) +
    theme(legend.key = element_rect(colour = NA)) +
    labs(linetype = ystrataname) +
    theme(plot.margin = unit(c(0, 1, .5, ifelse(m < 10, 1.5, 2.5)), "lines")) +
    ggtitle(main)
  
  if (pval) {
    sdiff <-
      survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <-
      pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
    pvaltxt <-
      ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
    p <-
      p + annotate("text", x = 0.03 * max(sfit$time), y = 0.1, label = pvaltxt)
  }
  
  ## Create a blank plot for place-holding
  ## .df <- data.frame()
  blank.pic <- ggplot(.df, aes(time, surv)) +
    geom_blank() +
    theme_bw() +
    theme(
      axis.text.x = element_blank(), axis.text.y = element_blank(),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.ticks = element_blank(), panel.grid.major = element_blank(),
      panel.border = element_blank()
    )
  if (table) {
    ## Create table graphic to include at-risk numbers
    risk.data <-
      data.frame(
        strata = summary(sfit, times = times, extend = TRUE)$strata,
        time = summary(sfit, times = times, extend = TRUE)$time,
        n.risk = summary(sfit, times = times, extend = TRUE)$n.risk
      )
    data.table <-
      ggplot(risk.data, aes(
        x = time, y = strata, label = format(n.risk, nsmall = 0)
      )) +
      #, color = strata)) +
      geom_text(size = 3.5) +
      theme_bw() +
      scale_y_discrete(breaks = as.character(levels(risk.data$strata)), 
                       labels = ystratalabs) +
      # scale_y_discrete(#format1ter = abbreviate,
      # breaks = 1:3,
      # labels = ystratalabs) +
      scale_x_continuous("Numbers at risk", limits = c(0, max(sfit$time))) +
      theme(
        axis.title.x = element_text(size = 10, vjust = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(face = "bold", hjust = 1)
      )
    data.table <- data.table + theme(legend.position = "none") +
      xlab(NULL) + ylab(NULL)
    data.table <- data.table +
      theme(plot.margin = unit(c(
        -1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5) - 0.28 * m
      ), "lines"))
    gridExtra::grid.arrange(
      p, blank.pic, data.table,
      clip = FALSE, nrow = 3, ncol = 1,
      heights = unit(c(2, .1, .25),c("null", "null", "null"))
    )
    a <- gridExtra::arrangeGrob(
      p, blank.pic, data.table, clip = FALSE,
      nrow = 3, ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null"))
    )
    return(a)
    }
  p
}