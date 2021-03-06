% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfmat.R
\name{removeColumns}
\alias{removeColumns}
\title{Remove columns}
\usage{
removeColumns(data.frame, columns, drop = FALSE)
}
\arguments{
\item{data.frame}{data.frame}

\item{columns}{names of columns to be removed}

\item{drop}{Logical, whether the matrix should be dropped to vector if only
one column is left}
}
\value{
data.frame with specified columns removed
}
\description{
Remove columns from a data.frame object
}
\details{
The function is equivalent to the subsetting operation with brackets. It
provides a tidy programming interface to manupulate data.frames.
}
\examples{

clubs <- data.frame(Points=c(21,23,28,24), Name=c("BVB", "FCB", "HSV",
"FCK"), games=c(12,11,11,12))
removeColumns(clubs,c("Name"))

}
\author{
Jitao David Zhang <jitao_david.zhang@roche.com>
}
