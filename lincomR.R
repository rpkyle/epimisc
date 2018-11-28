lincomR <- function(model, contrast, conflev = 0.95, digits = 3) {
	if ("geeglm" %in% class(model)) {
			modfamily <- model$family$family
			modlink <- model$family$link
			modcoef <- model$coeff
    		vcm <- model$geese$vbeta
		}
	else if ("geese" %in% class(model)) {
			modfamily <- model$model$variance
			modlink <- model$model$mean.link
			modcoef <- model$beta
			vcm = model$vbeta
		}
	else if ("gee" %in% class(model)) { 
			modfamily <- tolower(fm$model$varfun)
			modlink <- tolower(fm$model$link)
			modcoef <- model$coeff
			vcm <- model$robust.variance
		}
	else if ("glm" %in% class(model)) { 
			modfamily <- family(model)[[1]]
			modlink <- family(model)[[2]]			
			modcoef <- model$coeff
			vcm = vcov(model)
		}
	#else if ("glmmPQL" %in% class(model)) { 
	#		modfamily <- model$family$family
	#		modlink <- model$family$link
	#		modcoef <- model$coeff[[1]]
	#		vcm <- vcov(model)
	#	}
	else if ("coxph" %in% class(model)) { 
			modlink <- NA
			modcoef <- model$coeff
			vcm <- vcov(model)
		}
	else if (any(c("lmerMod", "glmerMod") %in% class(model))) {
			suppressPackageStartupMessages(library(lme4)) 
			modfamily <- family(model)[[1]]
			modlink <- family(model)[[2]]						
			modcoef <- model@beta
			vcm <- vcov(model)
		}
	else if ("glmmML" %in% class(model)) {
			modcoef <- model$coeff
			vcm <- model$variance[-nrow(model$variance),-ncol(model$variance)]
		}
	else if ("yagsResult" %in% class(model)) {
			modfamily <- model@family[[1]]
			modlink <- model@family[[2]]
			modcoef <- model@coefficients
			vcm <- model@robust.parmvar
		}
	else stop("Invalid model object type: lincomR currently supports model objects of class coxph, gee, geeglm, geese, glm, glmmML, mer, yags.")

	n = length(modcoef)

	if (n != length(contrast)) stop ("Contrast vector length not equal to number of model coefficients!")

	covar <- matrix(0,n,n)
	var <- matrix(0,n,n)

	for (i in 1:n) 

		logest <- sum(modcoef * contrast)
		expest <- exp(logest)

		for (i in 1:n) {
			for (j in 1:n) {
				if (i==j) 
					var[i,j] <- contrast[i]^2 * vcm[i,i] 
				else 
					covar[i,j] <- contrast[i] * contrast[j] * vcm[i,j]
			}
		}

	se <- sqrt( sum(var, na.rm=T) + sum(covar, na.rm=T) )
	se.expest <- exp(logest)*se

	zscore <- abs( qnorm( (1 - conflev)/2) )
	CI <- exp( logest + c(-zscore,zscore)*se)

	Z <- logest/se

	pvalue <- 2*pnorm(-abs(Z))
	
	if("yagsResult" %in% class(model))
		parmnames <- model@varnames
	else {
		parmnames <- names(modcoef)
	}
	
	parmvec <- contrast[contrast != 0]
	plusvec <- c(rep(" +", times = length(parmvec) - 1), "")

	cat("Linear combination of parameters:\n\n")
	cat(" ", paste0("(", parmvec, " * ", parmnames[1:n*as.numeric(contrast != 0)], 
		")", plusvec), "\n\n")

	if ("glmmML" %in% class(model)) {
		x <- data.frame(logest, se, Z, pvalue, expest, CI[1], CI[2])
		colnames(x) <- c("Estimate", " Std. Error", " Wald z", " P>|z|", " Exp(Est)", " Lower CL", " Upper CL")

		print(format(x, digits = digits), row.names = FALSE)		
		cat("\nGLM with random intercept fit using glmmML\n")
		cat("Note: lincomR assuming that link function is logit.")
	}
	else if ("coxph" %in% class(model)) {
		x <- data.frame(logest, se, Z, pvalue, expest, CI[1], CI[2])
		colnames(x) <- c("Estimate", " Std. Error", " Wald z", " P>|z|", " Exp(Est)", " Lower CL", " Upper CL")

		print(format(x, digits = digits), row.names = FALSE)		
		cat("\nCox proportional hazards model fit using coxph()\n")
		cat("Tie handling method:", model$method, "\n\n")
	}
	else if (any(c("geese", "glmmPQL", "yagsResult") %in% class(model))) {
		if(modlink %in% (c("logit", "log"))) {
			x <- data.frame(logest, se, Z, pvalue, expest, CI[1], CI[2])
			colnames(x) <- c("Estimate", " Std. Error", " Wald z", " P>|z|", " Exp(Est)", " Lower CL", " Upper CL")
		} else {
			CI <- logest + c(-zscore,zscore)*se
			x <- data.frame(logest, se, CI[1], CI[2], Z, pvalue)
			colnames(x) <- c("Estimate", " Std. Error", " Lower CL", " Upper CL", " Wald z", " P>|z|")
		}

		print(format(x, digits = digits), row.names = FALSE)
		
		cat("\nFamily:", modfamily, "\nLink function:", modlink, "\n\n")
	}
	else if (!("geese" %in% class(model))) {
		if(modlink %in% c("logit", "log")) {
			x <- data.frame(logest, se, Z, pvalue, expest, CI[1], CI[2])
			colnames(x) <- c("Estimate", " Std. Error", " Wald z", " P>|z|", " Exp(Est)", " Lower CL", " Upper CL")
		} else {
			CI <- logest + c(-zscore,zscore)*se
			x <- data.frame(logest, se, CI[1], CI[2], Z, pvalue)
			colnames(x) <- c("Estimate", " Std. Error", " Lower CL", " Upper CL", " Wald z", " P>|z|")
		}

		print(format(x, digits = digits), row.names = FALSE)		
		print(family(model))
	}

	if ("coxph" %in% class(model) | any(c("log", "logit") %in% modlink)) {
		result <- list(conflev = conflev, est = logest, se.est = se, wald.z = Z, pvalue = pvalue, exp.est = expest, se.exp.est = se.expest, lower.bound = CI[1], upper.bound = CI[2])
	} else {
		result <- list(conflev = conflev, est = logest, se.est = se, wald.z = Z, pvalue = pvalue, lower.bound = CI[1], upper.bound = CI[2])
	}
	invisible(result)
}