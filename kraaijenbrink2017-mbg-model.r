####  Dynamic mass-balance gradient glacier model with lumped mass redistribution component
####
####  Author:
####  Philip Kraaijenbrink (p.d.a.kraaijenbrink@uu.nl)
####
####  Reference:
####  Kraaijenbrink PDA, Bierkens MFP, Lutz AF and Immerzeel WW (2017).
####  Impact of a global temperature rise of 1.5 degrees Celsius on Asia’s glaciers. Nature.
####  http://doi.org/10.1038/nature23878
####
####  Required input data available for download at:
####  https://www.mountainhydrology.org/data-nature-2017
####
####


# Init =================================================================================================================

# clear the environment
rm(list=ls())

# load required R-packages
library(raster)
library(stringr)
library(optimx)
library(truncnorm)
library(parallel)
library(snowfall)
library(caTools)



# Settings ==============================================================================================================

# set workdir
setwd("/kraaijenbrink-mbg-model")

# input dirs
rdsroot  <- './input-data/rds-data'              # folder with RDS files
rastroot <- './input-data/raster-data'           # folder with RGIId-named subfolders with raster data

# output dirs
outroot <- "./output-data/mbg-model-output"
logroot <- "./output-data/mbg-model-logs"
dir.create(outroot,showWarnings=F, recursive=T)
dir.create(logroot,showWarnings=F, recursive=T)

# run multicore
runclustered  <- F
keepcoresfree <- 2



# Load and prep general data ===========================================================================================

# load centroid shape with precip, degree days, observed mb etc.
rgi.point     <- readRDS(file.path(rdsroot,'glacier-data.rds'))
rgi.ids       <- rgi.point$RGIId

# read lists with climate forcing (dTs and dPs for all models for all glaciers for all 5-year steps)
glac.pr.delta <- readRDS(file.path(rdsroot,"dP_factors_2006-2100.rds"))
glac.t.delta  <- readRDS(file.path(rdsroot,"dT_degrees_2006-2100.rds"))

# read the sampled mean ostrem curve
ostrem        <- readRDS(file.path(rdsroot,'ostrem_meancurve.rds'))

# make table with GCM run info
refyear       <- 2005
gcminfo       <- data.frame(model=gsub("_rcp.*","",gsub("dP_","",names(glac.pr.delta))),
                            rcp=str_extract(names(glac.pr.delta),("rcp.?.?")))




# Model function definition ================================================================================================

# debug glacier
counter <- 2594


# per-glacier function to call later in parallel scheme
mainCalcFunc <- function(counter){

  # glacier data of current iteration
  currglac <- rgi.point[counter,]@data

  # print progress to log file
  logtext <- sprintf('Processing iteration %05d, glacier %s', counter, currglac$RGIId)
  if (runclustered){flog.info(logtext)}else{print(logtext)}


  # Input parameter frame ----------------------------------------------------------------------------------------------------

  # initialize parameter data frame
  dfrows              <- 1+length(glac.pr.delta)
  gcmind              <- 2:(length(glac.pr.delta)+1)
  settings            <- data.frame(label=rep(NA,dfrows))
  settings$label      <- c('Current',  gsub('^dP_','',names(glac.pr.delta)))
  settings$gradient   <- 0
  settings$pr.ref     <- currglac$precip
  settings$dd.ref     <- currglac$dd_mean
  settings$ddf.cln    <- 7.0
  settings$debthick   <- 0.5
  settings$obs.mb.rnd <- 0

  # get indices of the projection ensembles
  cur.ind <-1
  r26.ind <- grep('rcp26',settings$label)
  r45.ind <- grep('rcp45',settings$label)
  r60.ind <- grep('rcp60',settings$label)
  r85.ind <- grep('rcp85',settings$label)
  det.ind <- c(cur.ind,r26.ind,r45.ind,r60.ind,r85.ind)


  # Monte Carlo parameters set initialisation -----------
  if (T){

    # initiate the random sample for the monte carlo
    modsamps       <- 10
    sample.n       <- length(glac.pr.delta)*modsamps
    ddf.cln.rand   <- runif(sample.n,0,1)
    ddf.deb.rand   <- runif(sample.n,0,1)
    dd.ref.rand    <- runif(sample.n,0,1)
    dP.rand        <- runif(sample.n,0,1)
    obs.mb.rand    <- runif(sample.n,0,1)
    debthick.rand  <- runif(sample.n,0,1)

    # DDF clean
    m              <- settings$ddf.cln[1]
    s              <- 2
    ddf.cln.set    <- qtruncnorm(ddf.deb.rand,a=0,b=Inf,mean=m,sd=s)

    # Degree days
    m              <- settings$dd.ref[1]
    s              <- currglac$dd_sd
    dd.ref.set     <- qnorm(dd.ref.rand, mean=m, sd=s)

    # delta precip
    dP.set         <- qtruncnorm(dP.rand, a=0,b=2.0, mean=1, sd=0.33)
    dP.set         <- ifelse(dP.set<1.0, 1/(2-dP.set),dP.set)

    # massbalance offset
    m              <- 0
    s              <- currglac$mb_sd
    obs.mb.set     <- qnorm(obs.mb.rand, mean=m, sd=s)
    # obs.mb.set   <- rep(settings$obs.mb[1],sample.n)

    # debris thickness set (for 95th percentile LST)
    m              <- settings$debthick[1]
    s              <- 0.2
    debthick.set   <- qtruncnorm(debthick.rand, a=0.01, b=5.0, mean=m, sd=s)

    # add monte carlo to settings frame
    settingssize                <- nrow(settings)
    settings                    <- settings[c(1:settingssize,rep(1,sample.n)),]
    mcind                       <- (1:sample.n)+settingssize
    gcmsamp                     <- rep(1:length(glac.pr.delta),each=sample.n/length(glac.pr.delta))
    settings$label[mcind]       <- sprintf('MC%04d %s', 1:sample.n, gsub('^dP_','',names(glac.pr.delta))[gcmsamp])
    settings$gradient[mcind]    <- 0
    settings$ddf.cln[mcind]     <- ddf.cln.set
    settings$debthick[mcind]    <- debthick.set
    settings$dd.ref[mcind]      <- dd.ref.set
    settings$pr.ref[mcind]      <- dP.set * settings$pr.ref[1]
    settings$obs.mb.rnd[mcind]  <- obs.mb.set

    # fix row names
    row.names(settings) <- 1:nrow(settings)

  }else{
    sample.n = 0
    mcind = NULL
  }

  # derivatives of settings
  settings$obs.mb   <- currglac$mb_mean + settings$obs.mb.rnd
  settings$acc.max  <- settings$pr.ref * 0.001
  settings$abl.max  <- -settings$dd.ref * settings$ddf.cln * 0.001
  settings$abl.max  <- ifelse(settings$abl.max > -1.0, -1.0, settings$abl.max)   # force ablation to be at least -1 m w.e. a-1



  # Read and process raster inputs --------------------------------------------------------------------------------------------

  dem.rast      <- raster(file.path(rastroot, currglac$RGIId,'srtm-elevation.tif'))
  slope.rast    <- raster(file.path(rastroot,currglac$RGIId,'srtm-slope.tif'))
  ice.thickness <- crop(raster(file.path(rastroot,currglac$RGIId,'ice-thickness.tif')), dem.rast)
  debris        <- raster(file.path(rastroot,currglac$RGIId,'classification.tif'))
  tir           <- raster(file.path(rastroot,currglac$RGIId,'ls8-composite-tsurf.tif'))
  cell.area     <- prod(res(dem.rast))

  # apply slope threshold to debris classification
  debslopethresh                    <- 24
  debris[slope.rast > debslopethresh & debris==2] <- 1
  missing.debris.pixels             <- (is.na(debris) == !is.na(dem.rast))
  debris[missing.debris.pixels]     <- 0
  ponds                             <- debris==3
  ice                               <- debris==1
  debris                            <- debris==2

  # subset TIR raster to only the debris covered part
  tir[debris==0]    <- NA
  tir               <- tir-cellStats(tir,min)

  # get dem stats
  dem.vals          <- na.omit(as.vector(dem.rast))
  elev.min          <- min(dem.vals,na.rm=T)
  elev.median       <- median(dem.vals,na.rm=T)
  elev.max          <- max(dem.vals,na.rm=T)
  dem.rast.px       <- length(dem.vals)

  # get stepsize and number of bands according to glacier size
  # limit bands to predefined range
  band_limits      <- c(6,25)
  scaler            <- 250
  elev.stepsize     <- (elev.max-elev.min) / round(dem.rast.px / scaler)
  elev.steps        <- seq(elev.min,elev.max,elev.stepsize)
  no.of.bands     <- length(elev.steps)-1
  if (no.of.bands<band_limits[1]){
    elev.stepsize   <- (elev.max-elev.min) / band_limits[1]
    elev.steps      <- seq(elev.min,elev.max,elev.stepsize)
    no.of.bands   <- length(elev.steps)-1
  }
  if (no.of.bands>band_limits[2]){
    elev.stepsize   <- (elev.max-elev.min) / band_limits[2]
    elev.steps      <- seq(elev.min,elev.max,elev.stepsize)
    no.of.bands   <- length(elev.steps)-1
  }

  # get the elevation bands from dem spatially
  # check if there are any bands without pixels or with less than 3 pixels
  # recalculate with less bands if so
  emptybands=T
  while (emptybands){
    dem.rast.cls <-cut(dem.rast,c(-Inf,elev.steps[-c(1,length(elev.steps))],Inf))
    emptybands <- (nrow(na.omit(freq(dem.rast.cls))) != no.of.bands)  |  (length(which(freq(dem.rast.cls)[,2] < 3)) > 0)
    if (emptybands){
      no.of.bands  <- no.of.bands-1
      elev.stepsize  <- (elev.max-elev.min) / no.of.bands
      elev.steps     <- seq(elev.min,elev.max,elev.stepsize)
    }
  }

  # get ice thickness and debris stats for each band
  compareframe <- data.frame(band     = 1:no.of.bands,
                             thickness = zonal(ice.thickness,dem.rast.cls)[,2],
                             debris    = zonal(as.integer(debris),dem.rast.cls)[,2],
                             pond      = zonal(as.integer(ponds),dem.rast.cls)[,2])
  compareframe$ice <- 1-(compareframe$debris+compareframe$pond)

  # fix rare zero-thickness case (e.g. for a band that only consists of a 2 cell wide corridor)
  if (sum(compareframe$thickness==0)){
    minthick <- min(compareframe$thickness[-which(compareframe$thickness==0)])
    if (minthick < 1){  # assign the thickness of the band with the min thickness to the 0 thickness band
      compareframe$thickness[which(compareframe$thickness==0)] <- minthick
    }else{              # if minimum thickness is larger than 1 m, assign 1 m as thickness
      compareframe$thickness[which(compareframe$thickness==0)] <- 1
    }
  }



  # get debris thickness and DDF from TIR -----------------------------------------------------------------------------------------
  # using exponential relation between TIR and DCT
  # Monte Carlo thicknesses are set to be equal to the TIR 95th percentile
  tir.frame  <- na.omit(cbind(unlist(as.data.frame(tir)),unlist(as.data.frame(dem.rast.cls))))
  if (nrow(tir.frame)>0){
    aggrcls     <- sort(unique(tir.frame[,2]))
    tir.95p     <- quantile(tir, probs=0.95);
    if (tir.95p < (0.25*cellStats(tir,max))){
      tir.95p = cellStats(tir,max)
    }
    if (tir.95p ==0){
      tir.95p = 5
    }
    dct.frame   <- tir.frame[,rep(1,nrow(settings))]
    exp.scaler  <- matrix(log(settings$debthick*100)*(1/tir.95p),ncol=nrow(settings), nrow=nrow(tir.frame), byrow=T)
    dct.frame   <- exp(exp.scaler * dct.frame)
    dct.frame   <- ifelse(dct.frame > 500, 500, dct.frame)     # limit debris thickness to 5 m
    dct.cls     <- aggregate(dct.frame, by=list(tir.frame[,2]), FUN=function(x) mean(x))[,-1]
    ddf.scl     <- matrix(approx(ostrem,xout=dct.frame)$y, nrow=nrow(tir.frame))
    ddf.scl.cls <- aggregate(ddf.scl, by=list(tir.frame[,2]), FUN=function(x) mean(x))[,-1]
    ddf.deb     <- ddf.scl.cls  * rep(settings$ddf.cln,each=nrow(ddf.scl.cls))
  }



  # calculate effective DDF per band -----------------------------------------------------------------------------------------

  # DDFtot = (DDFice*Aice + DDFdeb*Adeb + DDFpond*Apond)  / Atot
  ddf.pond <- approx(ostrem, xout=50)$y*10 * settings$ddf.cln   # ponds have 10 times the melt of thick debris, i.e. 50 cm.
  ddf.pond <- matrix(ddf.pond, nrow=no.of.bands, ncol=nrow(settings), byrow=T)
  ddf.ice  <- matrix(settings$ddf.cln, nrow=no.of.bands, ncol=nrow(settings), byrow=T)

  ddf.ice.eff  <- ddf.ice*compareframe$ice
  ddf.pond.eff <- ddf.pond*compareframe$pond
  ddf.deb.eff  <- ddf.ice.eff*0
  if (nrow(tir.frame)>0){  # only fill if there is debris, otherwise keep at 0 DDF
    ddf.deb.eff[aggrcls,] <- as.matrix(ddf.deb*compareframe$debris[aggrcls])
  }

  # get effective DDF and resulting MB reduction factor
  ddf.effective <- ddf.ice.eff + ddf.deb.eff + ddf.pond.eff
  red.factor    <- 1 - (ddf.effective / ddf.ice)




  # Create Band MetaData Frame -------------------------------------------------------------------------------------------------
  bmdf <- data.frame(band=1:no.of.bands)

  # fill frame values
  bmdf$step.min  <- elev.steps[-length(elev.steps)]
  bmdf$step.max  <- elev.steps[-1]
  bmdf$step.cntr <- (bmdf$step.min+bmdf$step.max)/2
  if (!emptybands){
    bmdf$cells   <- freq(dem.rast.cls)[1:no.of.bands,2]
  }else{   # failsafe if there are bands without any pixels
    tmp                          <- merge(compareframe$band,freq(dem.rast.cls),by=1,all=T)[1:no.of.bands,]
    tmp[which(is.na(tmp[,2])),2] <- 1
    bmdf$cells                   <- tmp[,2]
  }
  bmdf$area     <- bmdf$cells * cell.area
  bmdf$depth    <- compareframe$thickness
  bmdf$volume   <- bmdf$area * bmdf$depth
  bmdf$debris   <- compareframe$debris
  bmdf$pond     <- compareframe$pond


  # initiate the 3D Dynamic Model Calculation Array, used for vectorized MC calculations
  dn                   <- list(1:no.of.bands,
                               c('area','volume','deb.reduct','mb','mb.debfac','mb.m3','gradient','ice.flux','net.flux'),
                               settings$label)
  dmca                 <- array(data=NA,dim=sapply(dn,length),dn)
  dmca[,'deb.reduct',] <- red.factor
  dmca[,'deb.reduct',][which(dmca[,'deb.reduct',] > 1)]  <- 1-dmca[,'deb.reduct',][which(dmca[,'deb.reduct',] > 1)]




  # Volume area relation sampling -------------------------------------------------------------------------------------------------
  # Determine the volume area relation per elevation bands

  # no of samples per band
  vas.samples <- 20

  # get raster values as numeric vector
  bin.df   <- as.vector(dem.rast.cls)
  vol.df   <- as.vector(ice.thickness * prod(res(ice.thickness)))

  # get volume values per band and discard NA
  vol.list <- lapply(1:no.of.bands, function(x) na.omit(vol.df[bin.df==x]))

  # sample area and volume for different ice thicknesses
  vasCalc <- function(x){
    x.range   <- range(x)
    samp.vols <- seq(x.range[1],x.range[2],diff(x.range)/vas.samples)
    out.vols  <- c(0,sapply(samp.vols,function(i) sum(x[x<=i],na.rm=T))[-1])
    out.areas <- c(0,sapply(samp.vols,function(i) length(which(x<=i))*cell.area)[-1])
    return(data.frame(thickness=samp.vols/cell.area,area=out.areas,vol=out.vols))
  }
  vasframes <- lapply(vol.list, vasCalc)

  # if bands do not have enough unique area-vol combinations to fit a spline
  # (sometimes for very small bands in the acc zone), use linear relation
  uniquevasrows <- sapply(vasframes, function(x) nrow(unique(x[,2:3],margin=2)))
  vasframes     <- lapply(1:length(vasframes), function(x){
    if(uniquevasrows[x] <= 6){
      vasframes[[x]] <- data.frame(thickness=seq(0,bmdf$depth[x],bmdf$depth[x]/vas.samples),
                                   area=seq(0,bmdf$area[x],bmdf$area[x]/vas.samples),
                                   vol=seq(0,bmdf$volume[x],bmdf$volume[x]/vas.samples))
    }
    uniquerows <- as.numeric(row.names(unique(vasframes[[x]][,2:3],margin=2)))
    return(vasframes[[x]][uniquerows,])})

  # fit spline through samples to get area predictor
  vas.splines <- lapply(vasframes, function(x) smooth.spline(x$vol,x$area, spar=0.05))




  # Determine the mass balance gradient --------------------------------------------------------------------------------------------

  # optimize mass balance gradient using observed mass balance for the current climate
  # optimizer function
  mb.gradient <- rep(NA,nrow(settings))
  for (g in which(settings$gradient==0)){
    abl.max.cur <- settings$abl.max[g]
    acc.max.cur <- settings$acc.max[g]
    obs.mb.cur  <- settings$obs.mb[g]
    deb.red.cur <- dmca[,'deb.reduct',g]
    gradopt <- function(mbg){
      mb          <- c(abl.max.cur, abl.max.cur + bmdf$band * mbg * elev.stepsize)[1:no.of.bands]
      mb          <- ifelse(mb > acc.max.cur, acc.max.cur, mb)
      mb.m3       <- ifelse(mb>0, mb*bmdf$area, (mb-mb*deb.red.cur)*bmdf$area)
      squareddiff <- ((sum(mb.m3)/sum(bmdf$area)-obs.mb.cur)^2) + 0.00001*mbg
      return (squareddiff)
    }

    # Fit gradient for current run only. Increase precip if ELA is below threshold, constrain correction upper limit
    if (g==1){
      mbfit                <- -1
      ela                  <- -1
      ela.thresh           <- quantile(dem.vals, probs=c(0.25))
      settings$pr.ref.orig <- settings$pr.ref
      while((tail(mbfit,1)<0) | (ela<ela.thresh[1])){
        curgrad <- optim(0.01, fn=gradopt, method='Brent', lower=-1, upper=1)$par
        mbfit <- c(abl.max.cur,abl.max.cur+elev.stepsize*bmdf$band*curgrad)[1:no.of.bands]
        if (tail(mbfit,1)<0){
          obs.mb.cur <- obs.mb.cur*0.9
        }
        elaFinder <- function(x){abs(settings$abl.max[1] + elev.stepsize*curgrad*x)}
        ela.step  <- optimize(elaFinder,c(-1e5,1e5))$minimum
        ela       <- bmdf$step.cntr[1]+ elev.stepsize*ela.step
        if (ela < ela.thresh[1]) {acc.max.cur <- acc.max.cur*1.1}  # increase precip iteratively with 10%
        if (acc.max.cur > 3)  {acc.max.cur <- 3; break}            # limit corrected precip to 3000 mm
      }
      # update settings and monte carlo members with optimized precip
      settings$pr.ref[c(cur.ind,gcmind)]    <- acc.max.cur * 1000
      if(sample.n>0){settings$pr.ref[mcind] <- dP.set * settings$pr.ref[1]}
      settings$acc.max                      <-  settings$pr.ref * 0.001
    }

    # fit mb gradient for all projections and members
    # lower the observed mass balance until there is at least one band in the acc zone (i.e. force positive gradient)
    mbfit <- -1
    while(tail(mbfit,1)<0){
      mb.gradient[g] <- optim(0.01, fn=gradopt, method='Brent', lower=-1, upper=1)$par
      mbfit <- c(abl.max.cur,abl.max.cur+elev.stepsize*bmdf$band*mb.gradient[g])[1:no.of.bands]
      if (tail(mbfit,1)<0){
        obs.mb.cur <- obs.mb.cur*0.9
      }
    }
    settings$obs.mb[g] <- obs.mb.cur
  }
  # set the gradient for settings that do not use their own optimized gradient
  mb.gradient[which(!settings$gradient==0)] <- mb.gradient[settings$gradient[which(!settings$gradient==0)]]
  # set the altered abl.max in case
  settings$obs.mb[which(!settings$gradient==0)] <- settings$obs.mb[settings$gradient[which(!settings$gradient==0)]]

  # get the current ELA from fitted mb gradient, i.e. the zero mass balance crossing
  elaFinder <- function(x){abs(settings$abl.max[1] + elev.stepsize*mb.gradient[1]*x)}
  ela.step  <- optimize(elaFinder,c(-1e5,1e5))$minimum
  ela       <- bmdf$step.cntr[1]+ elev.stepsize*ela.step




  # Get annual climate change forcing ----------------------------------------------------------------------------------------------

  # load the dp and dt for current glacier for all models
  dPs                 <- t(cbind(1,do.call(rbind,lapply(glac.pr.delta,function(x) x[counter,]))))
  dTs                 <- t(cbind(0,do.call(rbind,lapply(glac.t.delta,function(x) x[counter,]))))
  meanyears           <- c(2005,seq(2007.5,2097.5,5))

  # do a running mean low pass filter (30 years window?)
  window              <- 6
  dPs.rm              <- apply(dPs,2,function(x) runmean(x,k=window))
  dTs.rm              <- apply(dTs,2,function(x) runmean(x,k=window))

  # make sure the first year has no CC effect yet
  dPs.rm[1,]          <- 1
  dPs.rm[2,]          <- colMeans(dPs.rm[c(1,3),])
  dTs.rm[1,]          <- 0
  dTs.rm[2,]          <- colMeans(dTs.rm[c(1,3),])

  # fit a very slightly smoothing spline through the averaged data to obtain a predictor model
  all.years           <- min(meanyears):2100
  smoothing           <- 0.05
  dPs.rmsm.model      <- lapply(1:ncol(dPs.rm), function(x) smooth.spline(meanyears,dPs.rm[,x],spar=smoothing))
  dTs.rmsm.model      <- lapply(1:ncol(dTs.rm), function(x) smooth.spline(meanyears,dTs.rm[,x],spar=smoothing))

  # predict values for every year using the predictor and put them in a data frame
  dPs.rmsm            <- do.call(cbind,lapply(1:length(dPs.rmsm.model),function(x) predict(dPs.rmsm.model[[x]],all.years)$y))
  dTs.rmsm            <- do.call(cbind,lapply(1:length(dTs.rmsm.model),function(x) predict(dTs.rmsm.model[[x]],all.years)$y))

  # add yearly values for the non-changing curent conditions
  dPs.rmsm            <- as.data.frame(cbind(1,dPs.rmsm))
  dTs.rmsm            <- as.data.frame(cbind(0,dTs.rmsm))

  # repeat columns to match MC models in the settings frame
  dPs.rmsm            <- dPs.rmsm[,c(1, 2:(length(glac.pr.delta)+1), rep(2:(length(glac.pr.delta)+1), each=sample.n/length(glac.pr.delta)))]
  dTs.rmsm            <- dTs.rmsm[,c(1, 2:(length(glac.pr.delta)+1), rep(2:(length(glac.pr.delta)+1), each=sample.n/length(glac.pr.delta)))]

  # set row and column names
  row.names(dPs.rmsm) <- all.years
  names(dPs.rmsm)     <- settings$label
  row.names(dTs.rmsm) <- all.years
  names(dTs.rmsm)     <- settings$label




  # Run model dynamically with a yearly timestep ------------------------------------------------------------------------------------------------

  # years to run
  maxyears <- 95

  # initiate Transient Volume and Area Array
  tvaa <- array(NA,dim=c(maxyears+1, 2 ,nrow(settings)), dimnames=list(refyear:(refyear+maxyears), c('volume','area'),settings$label))

  # set initial volume to glabtop modelled volumes
  dmca[,'area',]   <- bmdf$area
  dmca[,'volume',] <- bmdf$volume

  # start dynamic calculations
  for (t in refyear:(refyear+maxyears)){

    # model timestep
    tm <- (t-refyear)

    # Yearly CC and derivs ---------
    # keep climate stable if model is run beyond 2100
    if (t <= max(as.numeric( row.names(dPs.rmsm)))){
      dp.t <- unlist(dPs.rmsm[tm+1,])
      dt.t <- unlist(dTs.rmsm[tm+1,])
    }else{
      dp.t <-  unlist(tail(dPs.rmsm,1))
      dt.t <-  unlist(tail(dTs.rmsm,1))
    }

    # get the ela for current timestep
    deladt.lims <- c(56.5, 202.5)  # limit to 90pct range of Shea and Immerzeel
    deladt      <- 16.46*(settings$pr.ref.orig-379.61)^0.32
    deladt      <- ifelse(deladt<deladt.lims[1] | is.na(deladt),deladt.lims[1],deladt)
    deladt      <- ifelse(deladt>deladt.lims[2] | is.na(deladt),deladt.lims[2],deladt)
    ela.t       <- dt.t * deladt
    ela.t.abs   <- ela+ela.t

    # Calculate accumulation and ablation ---------
    accmax.t            <- settings$pr.ref * 0.001 * dp.t
    ablmax.t            <- settings$abl.max - mb.gradient * ela.t
    dmca[,'mb',]        <- unlist(lapply(1:nrow(settings),function(x) c(ablmax.t[x],ablmax.t[x]+elev.stepsize*
                                                                         bmdf$band*mb.gradient[x])[1:no.of.bands]))
    dmca[,'mb',]        <- ifelse(dmca[,'mb',] > matrix(accmax.t,ncol=nrow(settings),nrow=no.of.bands,byrow=T),
                                 matrix(accmax.t,ncol=nrow(settings),nrow=no.of.bands,byrow=T), dmca[,'mb',])
    # increase neg. mass balance of the lower band if it's used to simulate advancing below that band,
    # i.e. growth beyond initial volume. It's done by linear scaling with mb.gradient.
    dmca[1,'mb',]       <- ifelse(dmca[1,'volume',]>bmdf$volume[1],
                                 mean(c(dmca[1,'mb',],
                                        dmca[1,'mb',] - (mb.gradient*elev.stepsize) * ((dmca[1,'volume',]-bmdf$volume[1]) / mean(bmdf$volume)))),
                                 dmca[1,'mb',])
    dmca[,'mb.debfac',] <- ifelse(dmca[,'mb',]>0, dmca[,'mb',],
                                 dmca[,'mb',] - (dmca[,'deb.reduct',]) * dmca[,'mb',])
    dmca[,'mb.m3',]     <- dmca[,'area',] * dmca[,'mb.debfac',]




    # Determine reology parameter -----------
    # Only determine parameter on 1st timestep using minimisation function.
    # Get the rheology parameter that would result in the modelled ice thicknesses,
    # given the glacier length and all boundary conditions.
    if (t==refyear){

      rheology <- rep(NA,nrow(settings))
      rheology.init  <- 1.2e-8      # (1.2e-8 used by Marshall 2011)
      for (g in which(settings$gradient==0)){

        acczone        <- c(which(dmca[,"mb.debfac",g]>=0))
        ablzone        <- c(which(dmca[,"mb.debfac",g]<0))

        # outflux required to keep the acc zone bands of equal size
        equiflux.acc   <- rev(cumsum(rev(dmca[acczone,'mb.m3',g])))

        # yearly influx for bottom band (as fraction of local MB)
        curmb   <- dmca[ablzone,'mb.m3',g]                        # net mass balance for ablation zone, i.e. flux + mass loss
        totmb   <- settings$obs.mb[g]*sum(bmdf$area)              # total negative mass balance in m3

        # get the complete equiflux
        if (length(ablzone)==1){
          equiflux.abl    <- 0
        }else if (length(ablzone)>=2){
          # get the required influx by attenuation/amplification of mb by the actual mb
          attmb           <- ((totmb / sum(curmb)) * curmb)[1:(length(ablzone)-1)]
          equiflux.abl    <- c(0, -cumsum(curmb[1:length(attmb)]-attmb))
        }

        # combine equifluxes
        equiflux        <- c(equiflux.abl,equiflux.acc)           # combine acc and abl fluxes
        names(equiflux) <- 1:no.of.bands
        equiflux[1]     <- equiflux[2]                            # set the flux of the bottom band the same as the band above
                                                                  # (just to be able to calculate a glacier gradient for the the bottom band)

        # optimize the rheology, given the estimated equiflux and the RGI glacier length
        rheologyOpt <- function(x){
          gradient     <- (1/(x * bmdf$area * bmdf$depth^5 * equiflux^-1)) ^ (1/3)        # required gradient given the input rheology
          glaclength   <- elev.stepsize/ gradient                                         # glacier length given the modelled gradient
          sqlendiff    <- (sum(glaclength)-currglac$Lmax)^2                               # squared diff between modelled and "real" length
          return(sqlendiff)
        }
        rheopt         <- optim(rheology.init, rheologyOpt, method='Brent', lower=0, upper=1)  # run optimizer function
        rheology[g]    <- rheopt$par

        # run the gradient function with optimized rheology parameter and
        # fill the helper frame with estimated gradient and length values
        dmca[,'gradient',g] <- (1/(rheology[g] * bmdf$area * bmdf$depth^5 * equiflux^-1)) ^ (1/3)
      }
      # set the rheology for remainder of parameter combs that do not have their own mb gradient
      rheology[which(!settings$gradient==0)] <- rheology[settings$gradient[which(!settings$gradient==0)]]
      # set the glacier gradients for remainder of parameter combs that do not have their own mb gradient
      dmca[,'gradient',which(!settings$gradient==0)] <- dmca[,'gradient',settings$gradient[which(!settings$gradient==0)]]
    }
    # put into correct matrix representation for the vectorized Monte Calro calculations
    rheologymat <- matrix(rep(rheology,no.of.bands), ncol=nrow(settings), byrow=T)


    # Glacier flux  -----------

    # calculate the glacier flux using the current band volume with half of the
    # mass balance of that year to account for intra-annual dynamics
    dmca[,'volume',]     <- dmca[,'volume',] + 0.5*dmca[,'mb.m3',]
    for (c in 1:no.of.bands){dmca[c,'area',] <- predict(vas.splines[[c]],x=dmca[c,'volume',])$y}

    # flux calculation
    H                    <- ifelse((dmca[,'area',]>0) & (dmca[,'volume',]>0),                # Get average height of band
                                  dmca[,'volume',]/dmca[,'area',], 0)
    dmca[,'ice.flux',]   <- rheologymat * dmca[,'area',] * H^5 * dmca[,'gradient',]^3        # Calculate the downflux using Qi = rheo * Ai * Hi^5 * grad(Hi)^3
    dmca[1,'ice.flux',]  <- 0                                                                # set flux of bottom band back to zero
    influx               <- rbind(dmca[-1,'ice.flux',],0)                                    # Get influx, i.e. outflux of band above
    dmca[,'net.flux',]   <- influx - dmca[,'ice.flux',]                                      # Get the net horizontal flux

    # put half of the mass balance back again
    dmca[,'volume',]     <- dmca[,'volume',] - 0.5*dmca[,'mb.m3',]


    # New ice volume  -------------
    bandbalance         <- dmca[,'net.flux',] + dmca[,'mb.m3',]

    # limit band balance to a percentage of original band volume
    # any positive or negative excess is distributed iteratively over ice-containing bands
    # code runs until all excess have been dealt with
    relaxpar <- 0.2                                                                 # threshold
    non.empty <- apply(dmca[,'volume',],2,function(x) x>0)                          # bands with ice present
    excesscols <- which(apply(abs(bandbalance) > relaxpar*bmdf$volume,2,sum)>0)    # projs that have bands with excess
    for (i in excesscols){                                                          # loop over projs that have excess
      relaxexc <- c(1,1)                                                            # init excess for while loop
      while (sum(abs(relaxexc))>0){
        relaxbal <- ifelse((abs(bandbalance[,i]) > relaxpar*bmdf$volume),          # relaxed volume
                           relaxpar*bmdf$volume*sign(bandbalance[,i]),
                           bandbalance[,i])
        relaxexc <- bandbalance[,i] - relaxbal                                     # excess volume
        if (sum(abs(relaxexc))>0){
          topind   <- tail(which(abs(relaxexc[])>0),1)                              # topmost band with excess
          topexc   <- relaxexc[topind]                                              # topmost excess volume
          if (topind>1){                                                            # if excess is not in bottom band:
            distinds <- non.empty[,i]
            distinds[topind:length(relaxexc)] <- F
            if (sum(distinds)==0 & sum(relaxexc)<0){                                # if there are only emtpy bands below and excess is negative ...
              break                                                                 # ... retain the excess for that band
            }else{
              bandbalance[topind,i]   <- bandbalance[topind,i]-topexc             # in any other case put the topmost excess in the band below
              bandbalance[topind-1,i] <- bandbalance[topind-1,i]+topexc
            }
          }else{
            break
          }
        }
      }
    }
    dmca[,'volume',]   <- dmca[,'volume',] + bandbalance                            # update volume using relaxed band balance
    dmca[,'volume',]   <- ifelse(dmca[,'volume',]<0, 0, dmca[,'volume',])

    # hard limit the maximum volume in a band (that is not the bottommost one) to 5 times it's initial volume
    # (almost never required, failsafe for very rare numerical instability)
    dmca[-1,'volume',] <- ifelse(dmca[-1,'volume',] > 5*bmdf$volume[-1], 5*bmdf$volume[-1], dmca[-1,'volume',])

    # calculate change in area using a prediction from the volume-area relation splines
    dmca[is.na(dmca)]  <- 0
    for (c in 1:no.of.bands){dmca[c,'area',] <- predict(vas.splines[[c]],x=dmca[c,'volume',])$y}
    dmca[,'area',] <- ifelse(dmca[,'volume',] <= 0, 0, dmca[,'area',])
    dmca[,'area',] <- ifelse(dmca[,'area',] > bmdf$area, bmdf$area, dmca[,'area',])
    dmca[,'area',] <- ifelse(dmca[,'area',] < 0, 0, dmca[,'area',])

    # if a glacier grows beyond intital volume in the lowest band, apply a volume-area relation that simulates the glacier's advance
    dmca[1,'area',]  <- ifelse(dmca[1,'volume',]>bmdf$volume[1],
                               bmdf$area[1] + (dmca[1,'volume',]-bmdf$volume[1]) * mean(bmdf$area)/mean(bmdf$volume),
                               dmca[1,'area',])
    # if a glacier grows beyond intital volume in other bands, use linear scaling
    dmca[-1,'area',] <- ifelse(dmca[-1,'volume',]>bmdf$volume[-1],
                               dmca[-1,'area',] + dmca[-1,'area',]*((dmca[-1,'volume',] / bmdf$volume[-1])-1)*0.33,
                               dmca[-1,'area',])

    # add this timesteps data to the yearframe
    tvaa[tm+1,'volume',] <- colSums((dmca[,'volume',]))
    tvaa[tm+1,'area',]   <- colSums((dmca[,'area',]))

  } # end of the dynamic loop



  # Determine debris and ELA stats -------
  ice.thickness*cell.area -> icevolume.total -> icevolume.debris -> icevolume.clean
  icevolume.debris[ice!=0] <- NA
  icevolume.clean[ice!=1]  <- NA
  debris.area              <- length(which(as.matrix(ice==0)))*cell.area
  clean.area               <- length(which(as.matrix(ice==1)))*cell.area
  total.area               <- debris.area + clean.area
  debris.volume            <- sum(as.matrix(icevolume.debris),na.rm=T)
  clean.volume             <- sum(as.matrix(icevolume.clean),na.rm=T)
  total.volume             <- sum(as.matrix(icevolume.total),na.rm=T)
  debrisratio              <- as.data.frame(matrix(c(total.area,debris.area,clean.area,
                                                   total.volume,debris.volume,clean.volume),nrow=1))
  names(debrisratio)       <- c('total.area','debris.area','clean.area',
                                'total.volume','debris.volume','clean.volume')

  # debris volume below and above ELA
  below.ela <- dem.rast < ela
  ela.stat  <- data.frame(tot.area = c(cellStats(below.ela,sum)*cell.area, cellStats(!below.ela,sum)*cell.area),
                          deb.area = c(cellStats(below.ela*debris,sum)*cell.area,cellStats((!below.ela)*debris,sum)*cell.area),
                          tot.vol  = c(cellStats(below.ela*icevolume.total,sum),cellStats((!below.ela)*icevolume.total,sum)),
                          deb.vol  = c(cellStats(below.ela*icevolume.debris,sum),cellStats((!below.ela)*icevolume.debris,sum)))
  row.names(ela.stat) <- c("Below ELA", "Above ELA")




  # Output setup -----------------------------------------------------------------------------------

  # update settings frame for output
  settings$mb.grad <- mb.gradient

  # add current determ stats to the band metaframe
  bmdf$grad.cur <- dmca[,'gradient',cur.ind]
  bmdf$len.cur  <- elev.stepsize / bmdf$grad.cur

  # prepare metadata list
  infolist            <- list(NULL)
  infolist$rgi        <- as.character(currglac$RGIId)
  infolist$rasterpix  <- dem.rast.px
  infolist$no.bands   <- no.of.bands
  infolist$indices    <- list(sample.n=sample.n, cur.ind=cur.ind, gcmind=gcmind,
                           r15.ind=r26.ind, r45.ind=r45.ind, r60.ind=r60.ind, r85.ind=r85.ind,
                           mcind=mcind, gradfit.ind=which(settings$gradient==0))

  # write list of outputs to disk
  mbdat           <- list(meta <- NULL)
  mbdat$meta      <- infolist
  mbdat$settings  <- settings
  mbdat$bandvals  <- bmdf
  mbdat$debris    <- debrisratio
  mbdat$yr        <- tvaa
  mbdat$ela       <- list(ela.volumes=ela.stat, ela.final=as.data.frame(ela.t.abs))
  mbdat$rheology  <- rheology
  dir.create(file.path(outroot,currglac$RGIId),showWarnings=F,recursive=T)
  outfile         <- file.path(outroot,currglac$RGIId,'model-output.rds')
  saveRDS(mbdat,file=outfile)

  # model function returns no variables
  return(NULL)

} # end of the model function ================================================================================================
#  ===========================================================================================================================






# Function call ===============================================================================================================

# glacier to run, i.e. the indices of 'rgi.point' to run the model for
glacs_to_run <- 1:length(rgi.point)
glacs_to_run <- 1234

# randomize glacier order to be able to make meaningfull plots early on
glacs_to_run <- as.integer(sample(as.character(glacs_to_run)))

# run the model
if (runclustered){ # run function clustered on multiple cores

  # number of workers
  nworkers <- detectCores()-keepcoresfree

  # run model function parallelized
  sfInit(parallel=T,cpu=nworkers, type="SOCK")                                        # initialize cluster
  loadlib <- lapply(c('maptools','raster','futile.logger','truncnorm','caTools'),     # load required packages in the cluster
                    function(x) sfLibrary(x, character.only=T))
  sfExportAll()                                                                       # transfer model variables to the workers
  sfClusterSetupRNG()                                                                 # setup random number generator

  # initiate logfiles
  loginit <- function(logfile) {
    library(futile.logger)
    flog.appender(appender.file(logfile))
    NULL
  }
  sfClusterApply(sprintf(file.path(logroot,'mgb-model-worker-%02d.log'), seq_len(nworkers)), loginit)

  # run function and get total process time
  tout <- system.time(out <- sfLapply(glacs_to_run, mainCalcFunc))

  # stop cluster and print process time
  sfStop()
  print(tout)


}else{ # run regularly on single core
  tout <- system.time(sapply(glacs_to_run, mainCalcFunc))
  print(tout)
}


# EOF

