## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(chronosphere)

## ----package, echo= TRUE, eval=FALSE------------------------------------------
#  library(chronosphere)

## ----dat, echo= TRUE----------------------------------------------------------
data(dems)
data(clim)

## ----dems, echo= TRUE---------------------------------------------------------
dems

## ----dems2, echo= TRUE--------------------------------------------------------
proxy(dems)

## ----clim, echo= TRUE---------------------------------------------------------
clim

## ----climprox, echo= TRUE-----------------------------------------------------
proxy(clim)

## ----constructor, echo= TRUE--------------------------------------------------
# a stack of rasters
stackOfLayers <- dems@stack
# an index object
ind <- 1:10
names(ind) <- letters[1:10]
# a RasterArray
nra  <- RasterArray(index=ind, stack=stackOfLayers)
nra

## ----comb, echo= TRUE---------------------------------------------------------
# one raster
r1 <- raster()
values(r1) <- 1
# same structure, different value
r2 <-raster()
values(r2) <- 2
comb <- combine(r1, r2)
comb

## ----rbind, echo= TRUE--------------------------------------------------------
# bind dems to itself
cbind(dems, dems)

## ----le, echo= TRUE-----------------------------------------------------------
length(dems)

## ----cle, echo= TRUE----------------------------------------------------------
nrow(clim)
ncol(clim)

## ----dime, echo= TRUE---------------------------------------------------------
dim(dems)

dim(clim)

## ----demname, echo= TRUE------------------------------------------------------
names(dems)

## ----climane, echo= TRUE------------------------------------------------------
colnames(clim)
rownames(clim)

## ----climrename, echo= TRUE---------------------------------------------------
clim2 <- clim
colnames(clim2) <- c("temp", "prec")
clim2

## ----climrenamedim, echo= TRUE------------------------------------------------
dimnames(clim2)[[1]] <- 1:10
clim2

## ----layernames, echo= TRUE---------------------------------------------------
layers(clim)

## ----cells, echo= TRUE--------------------------------------------------------
ncell(dems)
nvalues(dems)

## ----single, echo= TRUE, plot=TRUE, fig.height=5.5----------------------------
one <- dems[["dem_45"]]
mapplot(one, col="earth")

## ----dual, echo= TRUE---------------------------------------------------------
dems[[c(1,2)]]

## ----doublerep, echo= TRUE----------------------------------------------------
# copy
dem2 <- dems
dem2[["dem_0"]] <- dem2[["dem_5"]]

## ----doublecompare, echo= TRUE------------------------------------------------
# but these two are now the same
dem2[[1]]
dem2[[2]]

## ----dem30, echo= TRUE--------------------------------------------------------
dems["30"]

## ----demdrop, echo= TRUE------------------------------------------------------
dem30 <- dems["30", drop=FALSE]
class(dem30)

## ----beyond, echo= TRUE-------------------------------------------------------
dems[4:12]

## ----demna, echo= TRUE--------------------------------------------------------
demna <- dems
demna[3] <- NA

## ----cellsbu, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
# character is necessary, as the row named "2003" is necessary
one <- clim["2003", "bio1"]
mapplot(one)

## ----entirerow, echo= TRUE----------------------------------------------------
clim["2005", ]

## ----entirecol, echo= TRUE----------------------------------------------------
clim[,"bio12"]

## ----crop, echo= TRUE, plot=TRUE, fig.height=5.5------------------------------
# crop to Australia
ext <- extent(c(                
  xmin = 106.58,
  xmax = 157.82,
  ymin = -45.23,
  ymax = 1.14 
)) 

# cropping all DEMS (Australia drifted in)
au<- crop(dems, ext)

# select the first element
mapplot(au[1], col="earth")

## ----resam, echo= TRUE, plot=TRUE, fig.height=5.5-----------------------------
template <- raster(res=5)

# resample all DEMS
coarse <- resample(dems, template)

# plot an elemnt
mapplot(coarse["45"], col="earth")

## ----palette, fig.height=3, fig.width=8---------------------------------------
showPal()

## ----rlayer_plot, echo=TRUE, out.width='80%', fig.align="center"--------------
data(dems)
mapplot(dems[1])

#using an custom colour palette
mapplot(dems[1], col="earth", main="Using the earth palette")


## ----rarray_plot_def, echo=TRUE, out.width='80%', fig.align="center"----------
data(dems)
mapplot(dems)


## ----rarray_plot_leg, echo=TRUE, fig.align="center"---------------------------
data(dems)
mapplot(dems, col="ocean", legend=TRUE, legend.title="DEM", 
        plot.title=proxy(dems))

## ----na_layer, echo=TRUE, fig.align="center"----------------------------------
data(dems)

#create NA layer
dems[5] <- NA
is.na(dems)

mapplot(dems, col="coldhot")


## ----rarray_plot_ncol, echo=TRUE----------------------------------------------
# 4 columns
data(dems)
mapplot(dems, ncol=4)

## ----rarray_plot_multi, echo = TRUE, eval=FALSE-------------------------------
#  data(dems)
#  mapplot(dems, multi=TRUE)

## ----rarray_plot_nvar, echo=TRUE, fig.width=8, fig.height=4.5-----------------
data(clim)

mapplot(clim[1:2,], legend=TRUE)

## ----rarray_plot_rowlabs, echo=TRUE, fig.width=8, fig.height=4.5--------------
data(clim)

mapplot(clim[1:2,], col=c("coldhot", "wet"), 
        legend=TRUE, legend.title=c("Temperature", "Precipitation"))

