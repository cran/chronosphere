% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot,RasterArray,missing-method}
\alias{plot,RasterArray,missing-method}
\alias{plot,SpatialArray,missing-method}
\title{Shorthand for the plotting \code{\link[chronosphere:RasterArray-class]{RasterArray}} and \code{\link[chronosphere:SpatialArray-class]{SpatialArray}} objects}
\usage{
\S4method{plot}{RasterArray,missing}(x, y, ...)

\S4method{plot}{SpatialArray,missing}(x, y, ...)
}
\arguments{
\item{x}{A (\code{\link[chronosphere:RasterArray-class]{RasterArray}} or \code{\link[chronosphere:SpatialArray-class]{SpatialArray}}) Object to be plotted.}

\item{y}{Not implemented yet.}

\item{...}{Arguments passed to the \code{\link{mapplot}} function.}
}
\value{
None.
}
\description{
This \code{plot}, method executes the \code{\link{mapplot}} function on the \code{\link[chronosphere:RasterArray-class]{RasterArray}} or \code{\link[chronosphere:SpatialArray-class]{SpatialArray}}object.
}
\examples{
data(dems)
plot(dems)
}
