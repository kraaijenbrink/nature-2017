# aggregate the transient data to RGI regions
# author: philip kraaijenbrink

# load packages
library(maptools)
library(sp)
library(raster)
library(pbapply)
library(abind)

# clear environment
rm(list=ls())

# set workdir
setwd("/kraaijenbrink-mbg-model")

# set input root with model output
inroot <- "./output-data/mbg-model-output"

# root with RDS input files
rdsroot  <- './input-data/rds-data'

# get list with glaciers from the dirlist
dirlist <- list.files(inroot,full.names=T)
rgi.ids <- basename(dirlist)

# output root for aggregated output (one file per RGI region)
outroot <- "./output-data/mbg-model-output-aggregated"
dir.create(outroot,showWarnings=F,recursive=T)

# read vector data
rgi.point   <- readRDS(file.path(rdsroot,'glacier-data.rds'))
rgi.subreg  <- readRDS(file.path(rdsroot,'rgi-subregions.rds'))

# get region names and abbreviations from the subreg shapefile
regnames <- as.character(rgi.subreg$Primary_ID)
regnames <- gsub('C ','Central ',gsub('S ','South ',gsub('N ','North ',gsub('E ','East ',gsub('W ','West ',regnames)))))
regnames <- gsub(' \\(.*\\)$','',regnames)
regnames.abbr <- c('HIA','PAM',"WTS","ETS","WKL","EKL","QLS","INT","SET","HKU","KRK","WHL","CHL","EHL","HDS")

# get the points and ids separately per region
reg.points  <- lapply(1:length(rgi.subreg), function(x) rgi.point[rgi.subreg[x,],])
reg.ids     <- lapply(reg.points, function(x) as.character(x@data$RGIId))
reg.dirs    <- lapply(reg.ids, function(x) dirlist[which(rgi.ids %in% x)])
reg.count   <- sapply(reg.dirs, length)

# function to get the aggregated output
regionalAggregate <- function(i){

  print(sprintf('Processing region %02d: %s',i,rgi.subreg[i,]$Label))

  # intiate frames
  rds        <- readRDS(file.path(reg.dirs[[i]][1],'model-output.rds'))
  debrisstat <- rds$debris * 0
  meta       <- rds$meta
  settings   <- rds$settings
  elavol     <- rds$ela$ela.volumes * 0
  yr         <- rds$yr[,1:2,] * 0

  # iterate over glaciers in a region and read data into frames cumulatively
  looprange  <- 1:length(reg.dirs[[i]])
  pb         <- txtProgressBar(min=min(looprange),max=max(looprange),initial=0, style=3,width=50)
  for (j in looprange){
    setTxtProgressBar(pb, j, title = NULL, label = NULL)

    # update data
    rds        <- readRDS(file.path(reg.dirs[[i]][j],'model-output.rds'))
    debrisstat <- debrisstat + rds$debris
    yr         <- yr + rds$yr[,1:2,]
    elavol     <- elavol + rds$ela$ela.volumes
  }
  close(pb)

  # handle output frames
  yv         <- as.data.frame(yr[,1,])
  names(yv)  <- settings$label;
  ya         <- as.data.frame(yr[,2,])
  names(ya)  <- settings$label;
  debrisstat <- as.data.frame(debrisstat)

  # get percentage-wise change
  yvp  <- yv / debrisstat$total.volume * 100
  yap  <- ya / debrisstat$total.area * 100

  # write to output
  outlist <- list(volume=yv,vol.perc=yvp,
                  area=ya,area.perc=yap,
                  debris=debrisstat,
                  elavol=elavol,
                  meta=meta,settings=settings)
  saveRDS(outlist, file=file.path(outroot,sprintf('%03d-%s-aggregated.rds',i,rgi.subreg[i,]$Abbrev)))

  return(NULL) # return nothing
} # end of function

# function call for all regions
runinds <- 1:length(rgi.subreg)
tout <- system.time(out <- lapply(runinds, regionalAggregate))
print(tout)