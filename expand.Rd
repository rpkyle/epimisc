\name{expand}
\alias{expand}
\title{
expand
}
\description{
The \code{expand} function is similar to the function of the same name available in Stata. Given a dataset, and an expansion factor N, expand will output all rows in a data frame, repeated N times each.
}
\usage{
expand(dataset, size)
}
\arguments{
  \item{dataset}{
A data frame or matrix to expand.
}
  \item{size}{
Number of times each row should be output in the resulting data frame or matrix.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ryan Patrick Kyle <ryan.kyle@mail.mcgill.ca>
}
\examples{
## The function is currently defined as
function (dataset, size) 
{
    size <- abs(size)
    if (!(all(floor(size) == size, na.rm = TRUE) & length(size) == 
        1)) {
        stop("argument size must be an integer of length 1!")
    }
    return(dataset[rep(1:nrow(dataset), each = size), ])
  }
}
