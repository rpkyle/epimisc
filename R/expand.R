expand <-
function(dataset, size) {
	size <- abs(size)
	if(!(all(floor(size) == size, na.rm = TRUE) & length(size) == 1)) {
		stop("argument size must be an integer of length 1!")
	} # verify value is an integer of length one
	return(dataset[rep(1:nrow(dataset), each=size),])
}
