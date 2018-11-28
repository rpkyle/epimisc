cochran_qtest <-
function(estimate1, lcl1, ucl1, estimate2, lcl2, ucl2) {
	var1 <- (((log(ucl1)-log(estimate1))/qnorm(0.025)))^2
	var2 <- (((log(ucl2)-log(estimate2))/qnorm(0.025)))^2	
	num_pool <- (log(estimate1)/var1) + (log(estimate2)/var2)
	denom_pool <- (1/var1) + (1/var2)
	pooled <- num_pool/denom_pool
	q <- (((log(estimate1)-pooled)^2)/var1) + (((log(estimate2)-pooled)^2)/var2)
	return(list(pooled.hr = round(exp(pooled), 2), cochran.q = round(q, 2), P = round(1-pchisq(q, 1), 4)))
}
