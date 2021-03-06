% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inter.R
\name{matchtime}
\alias{matchtime}
\alias{matchtime,numeric-method}
\alias{matchtime,character-method}
\alias{matchtime,RasterArray-method}
\title{Match the dates of a time-dependent variable with a predefined vector}
\usage{
matchtime(x, y, ...)

\S4method{matchtime}{numeric}(x, y, index = FALSE, ...)

\S4method{matchtime}{character}(x, y, index = FALSE, ...)

\S4method{matchtime}{RasterArray}(x, y, index = FALSE, time = 1, ...)
}
\arguments{
\item{x}{Object to be reordered to match \code{y}.}

\item{y}{(\code{numeric}) The vector of dates (numeric values) to order to.}

\item{...}{Additional arguments passed to class-specific methods.}

\item{index}{(\code{logical}) If this argument is \code{TRUE}, only the indices will be returned that refer to the new order, rather than the reordered \code{x} variable.}

\item{time}{\code{numeric}. Single value referring to that dimension of \code{x} where the time-coding names are (\code{time=1} is the default for RasterArrays in \code{chronosphere}).}
}
\value{
An object of the class as \code{x} or a \code{numeric} vector.
}
\description{
The function takes a variable \code{x} (e.g. a vector or a \code{\link[chronosphere:RasterArray-class]{RasterArray}} object), and reorders it to best match the dates provided in a vector \code{y}.
}
\examples{
# original vector
orig <- 1:10
# target values
targ <- c(5.1,4.2, 3.4, 2.7, 2.3)
# how do the two series match the best?
matchtime(orig, targ)
}
