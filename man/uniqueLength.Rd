% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vector.R
\name{uniqueLength}
\alias{uniqueLength}
\alias{ulen}
\title{Length of unique elements in a vector}
\usage{
uniqueLength(x, incomparables = FALSE)
}
\arguments{
\item{x}{A vector}

\item{incomparables}{See \code{\link{unique}}}
}
\value{
An integer indicating the number of unique elements in the input
vector
}
\description{
Length of unique elements in a vector
}
\examples{

test.vec1 <- c("HSV", "FCB", "BVB", "HSV", "BVB")
uniqueLength(test.vec1)

test.vec2 <- c(1L, 2L, 3L, 5L, 3L, 4L, 2L, 1L, 5L)
ulen(test.vec2)

}
\seealso{
\code{unique}
}
\author{
Jitao David Zhang <jitao_david.zhang@roche.com>
}
