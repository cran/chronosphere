% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colours.R
\name{ramps}
\alias{ramps}
\alias{gradinv}
\alias{ocean}
\alias{terra}
\alias{ipccTemp}
\alias{ipccPrec}
\alias{wet}
\alias{ipccLine}
\alias{ipccRCP}
\alias{showPal}
\title{Colour gradient ramps}
\usage{
gradinv(n)

ocean(n)

terra(n)

ipccTemp(n, force = 11)

ipccPrec(n, force = 11)

wet(n)

ipccLine(n = 6)

ipccRCP(n = 4)

showPal(pal = "all")
}
\arguments{
\item{n}{(\code{numeric}) Number of different colors to generate from the palette}

\item{force}{(\code{logical}) Specify \code{pal} when multiple are available. More details to come.}

\item{pal}{(\code{character}) A palette name from the lists below}
}
\value{
A function producing a colour gradient ramp.
}
\description{
The object contains functions produced by the \code{\link[grDevices:colorRamp]{colorRampPalette}} function.
}
\details{
\code{showPal} can be used to display the available palettes. You can use \code{pal = "all"} or \code{pal=""} if you want to look at all the available palettes. 
You can also view single palettes individually. The following colour palettes are implemented:
\itemize{
\item \code{gradinv()}: inverse heatmap.
\item \code{ocean()}: ocean bathymetrical colours.
\item \code{terra()}: terrestrial topographical colours.
\item \code{ipccTemp()}: gradient from blue to red according to the official IPCC AR6 WG2 colour palette.
\item \code{ipccPrec()}: gradient from brown to green according to the official IPCC AR6 WG2  colour palette.
\item \code{wet()}: gradient from white to green to blue. 
\item \code{ipccLine()}: discrete colours for line graphs according to the official IPCC AR6 WG2  colour palette. 
\item \code{ipccRCP()}: discrete colours for climate scenarios according to the official IPCC AR6 WG2  colour palette. 
}
}
