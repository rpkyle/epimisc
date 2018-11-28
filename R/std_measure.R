std_measure <-
function(dataset, conflev=0.95, method, measure, user.weights = NA, digits = 3) {
  if (length(dim(dataset)) == 3) {
		  	n.strata <- dim(dataset)[3]
	  		dataset <- data.frame(y1=dataset[1,1,c(1:n.strata)], n1=colSums(dataset[1, c(1:2), c(1:n.strata)]), 
	  		y2=dataset[2,1,c(1:n.strata)], n3=colSums(dataset[2, c(1:2), c(1:n.strata)]))
	  } else if(!is.table(dataset)) stop("std.measure only accepts a table object as input")
	  else {
	  		n.strata <- 1
			dataset <- data.frame(y1=dataset[1,1], n1=sum(dataset[1, c(1:2)]), 
				y2=dataset[2,1], n3=sum(dataset[2, c(1:2)])) 
	}
  if (!method %in% c("MH", "internal", "external", "user.weights")) 
	  stop("Invalid method specified! Valid options are MH, internal, external, or user.weights.")
  if (!measure %in% c("risk.ratio", "risk.difference"))
	  stop("Invalid measure specified! Valid options are risk.ratio or risk.difference in this release.")	  
  if (method != "MH") {
	  if (!any(is.na(user.weights)) & method == "user.weights") {
		  if (length(user.weights) != nrow(dataset)) stop("Length of weights vector must equal number of strata")
	      W = user.weights; dataset <- cbind(dataset, W) 
	  }
	  else if (method == "internal") {
		  W = apply(dataset,1,function(ro) ro[2]); dataset <- cbind(dataset, W) 
		  if (any(!is.na(user.weights))) warning("User-specified weights are ignored if you do not also specify user as method.")
	  }
	  else if (method == "external") {
		  W = apply(dataset,1,function(ro) ro[4]); dataset <- cbind(dataset, W) 
		  if (any(!is.na(user.weights))) warning("User-specified weights are ignored if you do not also specify user as method.")
	  }
	  if (measure == "risk.difference") {
		  # standardized RD
		  delta = sum(apply(dataset,1,function(ro) ( ro[5]*( (ro[1]/ro[2]) - (ro[3]/ro[4]) ))))/sum(dataset[,5])
		  names(delta) <- "   Standardized Risk Difference"

		  # SE of standardized RD, and 95% CI for standardized RD
		  var.estimate <- (1/(sum(dataset[,5]))^2) * sum(apply(dataset,1,function(ro) ((ro[5]^2) * ( ((ro[1]*(ro[2]-ro[1]))/ro[2]^3) + ((ro[3]*(ro[4]-ro[3]))/ro[4]^3) ))))		  
		  se.delta = sqrt(var.estimate)		  
	      CI = delta + c(-1,1)*qnorm(1-(1-conflev)/2)*se.delta
	  } else if (measure == "risk.ratio") {
		  # standardized RR
		  delta = sum(apply(dataset,1,function(ro) (ro[5]*(ro[1]/ro[2]) ) )) / sum(apply(dataset,1,function(ro) (ro[5]*(ro[3]/ro[4]) ) ))
		  names(delta) <- "        Standardized Risk Ratio"
		  
		  # SE of standardized RR, and 95% CI for standardized RR
		  ln.var.estimate = (sum(apply(dataset, 1, function(ro) (ro[5]^2 * ro[1]) * ((ro[2] - ro[1])/(ro[2]^3)))) / (sum(apply(dataset, 1, function(ro) (ro[5]*(ro[1]/ro[2])) ))^2) ) +
	(sum(apply(dataset, 1, function(ro) (ro[5]^2 * ro[3]) * ((ro[4] - ro[3])/(ro[4]^3)))) / (sum(apply(dataset, 1, function(ro) (ro[5]*(ro[3]/ro[4])) ))^2) )		  
		  se.delta = exp(sqrt(ln.var.estimate))
		  CI=exp(log(delta) +c(-1,1)*qnorm(1-(1-conflev)/2)*(sqrt(ln.var.estimate)))
	  }
  }  
  else if (method == "MH" & measure == "risk.difference") {
	  num = sum(apply(dataset,1,function(ro) (ro[1]*ro[4] - ro[3]*ro[2])/(ro[2]+ro[4])))
	  if (any(is.na(user.weights))) {
		  W = apply(dataset,1,function(ro) ro[2]*ro[4]/(ro[2]+ro[4])) #Cochran weights  	
	  } else W = user.weights
	  dataset <- cbind(dataset, W)
	  
	  delta = num/sum(W)

      p1 = dataset[,1]/dataset[,2]
      p2 = dataset[,3]/dataset[,4]
      denom = apply(dataset,1,function(ro) ro[2]*ro[4]/(ro[2]+ro[4])) #Cochran weights
      var.MH.estimate = sum (denom^2 * (p1*(1-p1)/dataset[,2] + p2*(1-p2)/dataset[,4])) / sum(W)^2
	  se.delta = sqrt(var.MH.estimate)
	  CI=delta +c(-1,1)*qnorm(1-(1-conflev)/2)*sqrt(var.MH.estimate)
  } else if (method == "MH" & measure == "risk.ratio") {
	  num = sum(apply(dataset,1,function(ro) ((ro[1]*ro[4])/(ro[2]+ro[4]))))
	  if (any(is.na(user.weights))) {
		  W = apply(dataset,1,function(ro) ro[3]*ro[2]/(ro[2]+ro[4])) #Cochran weights
	  } else W = user.weights
	  dataset <- cbind(dataset, W) 
	  delta = num/sum(W)
					  
	  var.MH.num = sum(apply(dataset, 1, function(ro) ( ( (((ro[1]+ro[3]) * ro[2] * ro[4])/(ro[2]+ro[4])^2) - ((ro[1]*ro[3])/(ro[2]+ro[4])) ))))
	  var.MH.denom = sum(apply(dataset, 1, function(ro) ((ro[1]*ro[4])/(ro[2]+ro[4])))) * 
	  		sum(apply(dataset, 1, function(ro) ((ro[3]*ro[2])/(ro[2]+ro[4])))) 
	  ln.var.MH.estimate = var.MH.num/var.MH.denom
	  se.delta = exp(sqrt(ln.var.MH.estimate))
	  CI=exp(log(delta) +c(-1,1)*qnorm(1-(1-conflev)/2)*sqrt(ln.var.MH.estimate))	  
  	}

	cat("\nStandardization method:", method, "\n")
	cat("Specified alpha:", 1-conflev, "\n")
	
	if(measure == "risk.ratio") {
   	    # stratum-specific RRs
	    strdelta = (apply(dataset,1,function(ro) (ro[5]*(ro[1]/ro[2]) ) )) / (apply(dataset,1,function(ro) (ro[5]*(ro[3]/ro[4]) ) ))	  

		# SE of stratum-specific RRs, and 95% CI for stratum-specific RRs
		ln.var.strestimate = ((apply(dataset, 1, function(ro) (ro[5]^2 * ro[1]) * ((ro[2] - ro[1])/(ro[2]^3)))) / ((apply(dataset, 1, function(ro) (ro[5]*(ro[1]/ro[2])) ))^2) ) +
((apply(dataset, 1, function(ro) (ro[5]^2 * ro[3]) * ((ro[4] - ro[3])/(ro[4]^3)))) / ((apply(dataset, 1, function(ro) (ro[5]*(ro[3]/ro[4])) ))^2) )		  
		se.strdelta = exp(sqrt(ln.var.strestimate))		  		  
		strCI_lower = exp(log(strdelta) + -1 * qnorm(1-(1-conflev)/2)*(sqrt(ln.var.strestimate)))
		strCI_upper = exp(log(strdelta) + 1 * qnorm(1-(1-conflev)/2)*(sqrt(ln.var.strestimate)))		  

	    # crude RR
	    crudelta = (sum(dataset[1])/sum(dataset[2])) / (sum(dataset[3])/sum(dataset[4]))
		
		# SE of crude RR, and 95% CI for crude RR
		ln.var.cruestimate <- (sum(dataset[2]-dataset[1])/(sum(dataset[1])*sum(dataset[2]))) + (sum(dataset[4]-dataset[3])/(sum(dataset[3])*sum(dataset[4])))
		se.crudelta = exp(sqrt(ln.var.cruestimate))		  		  
		cruCI_lower = exp(log(crudelta) + -1 * qnorm(1-(1-conflev)/2)*(sqrt(ln.var.cruestimate)))
		cruCI_upper = exp(log(crudelta) + 1 * qnorm(1-(1-conflev)/2)*(sqrt(ln.var.cruestimate)))

		x1 <- data.frame(strdelta, se.strdelta, strCI_lower, strCI_upper, rownames(dataset), W)
		x2 <- data.frame(crudelta, se.crudelta, cruCI_lower, cruCI_upper)
		x3 <- data.frame(delta, se.delta, CI[1], CI[2])
		colnames(x1)[1] <- "        Stratum Risk Ratio"
		colnames(x2)[1] <- "          Crude Risk Ratio"
		
		if(method == "MH") colnames(x3)[1] <- "Mantel-Haenszel Risk Ratio"
			else colnames(x3)[1] <-          "   Standardized Risk Ratio"
	}
	else if(measure == "risk.difference") {		
		# stratum-specific RDs
		strdelta = (apply(dataset,1,function(ro) ( ro[5]*( (ro[1]/ro[2]) - (ro[3]/ro[4]) ))))/ (dataset[,5])
		
		# SE of stratum-specific RDs, and 95% CIs for stratum-specific RDs
		var.strestimate <- (1/((dataset[,5]))^2) * (apply(dataset,1,function(ro) ((ro[5]^2) * ( ((ro[1]*(ro[2]-ro[1]))/ro[2]^3) + ((ro[3]*(ro[4]-ro[3]))/ro[4]^3) ))))		  
		se.strdelta = sqrt(var.strestimate)
		strCI_lower = strdelta + -1 * qnorm(1-(1-conflev)/2)*(se.strdelta)
		strCI_upper = strdelta + 1 * qnorm(1-(1-conflev)/2)*(se.strdelta)		  
		
	    # crude RD
	    crudelta = (sum(dataset[1])/sum(dataset[2]))-(sum(dataset[3])/sum(dataset[4]))

		# SE of crude RD, and 95% CIs for crude RD
		var.cruestimate <- ((sum(dataset[1])*sum(dataset[2]-dataset[1]))/(sum(dataset[2])^3)) + ((sum(dataset[3])*sum(dataset[4]-dataset[3]))/(sum(dataset[4])^3))
		se.crudelta = sqrt(var.cruestimate)  
		cruCI_lower = crudelta + -1 * qnorm(1-(1-conflev)/2)*(se.crudelta)
		cruCI_upper = crudelta + 1 * qnorm(1-(1-conflev)/2)*(se.crudelta)		  

		x1 <- data.frame(strdelta, se.strdelta, strCI_lower, strCI_upper, rownames(dataset), W)
		x2 <- data.frame(crudelta, se.crudelta, cruCI_lower, cruCI_upper)
		x3 <- data.frame(delta, se.delta, CI[1], CI[2])
		colnames(x1)[1] <- "        Stratum Risk Difference"
		colnames(x2)[1] <- "          Crude Risk Difference"

		if(method == "MH") {
			colnames(x3)[1] <- "Mantel-Haenszel Risk Difference"
		} else colnames(x3)[1] <-          "   Standardized Risk Difference"
	}
	
	colnames(x1)[2:6] <- c(" Std. Error", " Lower CL", " Upper CL", "Stratum", "Weight")
	colnames(x2)[2:4] <- c(" Std. Error", " Lower CL", " Upper CL")
	colnames(x3)[2:4] <- c(" Std. Error", " Lower CL", " Upper CL")

	cat("\n")
	print(format(x1, digits = digits), row.names = FALSE)
	cat("\n")
	print(format(x2, digits = digits), row.names = FALSE)
	cat("\n")
	print(format(x3, digits = digits), row.names = FALSE)
	cat("\n")

	if(method == "MH" & measure == "risk.ratio") {
		pvalue <- (log(strdelta)-log(mhdelta))^2/ln.var.strestimate
		print(format(x3, digits = digits), row.names = FALSE)	  
	}
	result <- list(measure = measure, method = method, weights = W, conflev=conflev, std.estimate = unname(delta), se.std.estimate = se.delta, std.CI = CI)
	invisible(result)
}
