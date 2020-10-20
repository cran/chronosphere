## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(chronosphere)

## ----package, echo= TRUE, eval=FALSE------------------------------------------
#  library(chronosphere)

## ----dat, echo= TRUE----------------------------------------------------------
data(dems)

## ----dems, echo= TRUE---------------------------------------------------------
dems

## ----dems2, echo= TRUE--------------------------------------------------------
proxy(dems)

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

## ----coasts, echo= TRUE-------------------------------------------------------
data(coasts)
coasts

## ----coastprox, echo= TRUE----------------------------------------------------
proxy(coasts)

## ----combSp, echo= TRUE-------------------------------------------------------
# Individual SP objects
margin0 <- coasts[1, 1]
margin5 <- coasts[2, 1]
combination<- combine(margin0, margin5)
combination

## ----combre, echo= TRUE-------------------------------------------------------
coast0 <- coasts[1, 2]
coast5 <- coasts[1, 2]
# create SpatialStack
spStack<- stack(margin0, margin5, coast0, coast5)
# alternatively: accepts filenames too
# spStack <- SpatialStack(list(margin0, margin5, coast0, coast5))

## ----indcomb, echo= TRUE------------------------------------------------------
ind <- matrix(1:4, ncol=2)
colnames(ind) <- c("margin", "coast")
rownames(ind) <- c("0", "5")
spArray <- SpatialArray(index=ind, stack=spStack)
spArray


## ----le, echo= TRUE-----------------------------------------------------------
length(dems)

## ----cle, echo= TRUE----------------------------------------------------------
nrow(coasts)
ncol(coasts)

## ----dime, echo= TRUE---------------------------------------------------------
dim(dems)

dim(coasts)

## ----demname, echo= TRUE------------------------------------------------------
names(dems)

## ----coastname, echo= TRUE----------------------------------------------------
colnames(coasts)
rownames(coasts)

## ----coastrename, echo= TRUE--------------------------------------------------
coasts2 <- coasts
colnames(coasts2) <- c("m", "c")
coasts2

## ----caostrenamedim, echo= TRUE-----------------------------------------------
dimnames(coasts2)[[1]] <- paste(rownames(coasts), "Ma")
coasts2

## ----layernames, echo= TRUE---------------------------------------------------
layers(dems)

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
# character type is necessary as the rowname is "2003"
one <- coasts["15", "coast"]
mapplot(one)

## ----entirerow, echo= TRUE----------------------------------------------------
coasts["15", ]

## ----entirecol, echo= TRUE----------------------------------------------------
coasts[,"coast"]

## ----crop, echo= TRUE, plot=TRUE, fig.height=5.5------------------------------
# crop to Australia
ext <- extent(c(                
  xmin = 106.58,
  xmax = 157.82,
  ymin = -45.23,
  ymax = 1.14 
)) 

# cropping all DEMS 
au<- crop(dems, ext)

# select the first element
mapplot(au[1], col="earth")

## ----resam, echo= TRUE, plot=TRUE, fig.height=5.5-----------------------------
template <- raster(res=5)

# resample all DEMS
coarse <- resample(dems, template)

# plot an element
mapplot(coarse["45"], col="earth")

## ----sptran, echo= TRUE, plot=TRUE, fig.height=5.5----------------------------
# Cylindrial equal area projection
mollCoasts <- spTransform(coasts, CRS("+proj=moll"))

# plot edge of the map
edge<- spTransform(mapedge(), CRS("+proj=moll"))
plot(edge, col=ocean(8)[5])

# the entire thing
plot(mollCoasts["20", "margin"], col=ocean(8)[7], border=NA, add=T)
plot(mollCoasts["20", "coast"], col=terra(8)[5], add=T, border=NA)


## ----palette, fig.height=3, fig.width=8---------------------------------------
showPal()

## ----rlayer_plot, echo=TRUE, out.width='80%', fig.align="center"--------------
data(dems)
mapplot(dems[1])

#using a custom colour palette
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

mapplot(dems, col="ipccTemp")


## ----rarray_plot_ncol, echo=TRUE----------------------------------------------
# 4 columns
data(dems)
mapplot(dems, ncol=4)

## ----rarray_plot_multi, echo = TRUE, eval=FALSE-------------------------------
#  data(dems)
#  mapplot(dems, multi=TRUE)

