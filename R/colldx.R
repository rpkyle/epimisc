colldx <- function (mod, digits=3) {
  # Function written by Beau Bruce, based on collingenmodv9c.sas by Matthew Zack,
  # modified by Jim Singleton -- and colldiag from John Hendrickx's perturb package
  result <- NULL
  
  if ("geeglm" %in% class(mod)) {
    vcm <- mod$geese$vbeta
  }
  else if ("geese" %in% class(mod)) {
    vcm <- mod$vbeta
  }
  else if ("gee" %in% class(mod)) { 
    vcm <- mod$robust.variance
  }
  else if (any(c("glm", "lm", "coxph", "lmerMod", "glmerMod") %in% class(mod))) { 
    vcm <- vcov(mod)
  }
  else if ("glmmML" %in% class(mod)) {
    vcm <- mod$variance[-nrow(mod$variance),-ncol(mod$variance)]
  }
  else if ("yagsResult" %in% class(mod)) {
    vcm <- mod@robust.parmvar
  }
  else stop("Invalid model object type: colldx currently supports model objects of class coxph, gee, geeglm, geese, glm, glmmML, mer, yags.")
  
  ivcv <- solve(vcm)
  divcv <- diag(,ncol(vcm))*diag(ivcv)
  scale <- solve(sqrt(divcv))
  R <- scale %*% ivcv %*% scale
  svdR <- svd(R)
  val <- svdR$d
  vec <- svdR$v
  condindx <- sqrt(val[1]/val)
  Phi <- t((vec^2) %*% (diag(,length(val))*(1/val)))
  pi <- prop.table(Phi, 2)
  dim(condindx) <- c(length(condindx), 1)
  colnames(condindx) <- "cond.index"
  rownames(condindx) <- 1:nrow(condindx)
  colnames(pi) <- colnames(vcm)
  result$condindx <- condindx
  result$pi <- pi
  
  resdf <- as.data.frame(do.call(cbind, result))
  colnames(resdf)[1] <- "Condition Index"
  
  cat("Collinearity Diagnostics - Variance Decomposition Proportions\n\n")
  print(resdf, digits=digits, row.names = FALSE)
  invisible(result)
}
