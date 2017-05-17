##########################
#
# examine snow dynamics 
# for picker project
# preliminary agu analyses
# for 2016 fall mtg
#
# MML 11/23/16
##########################

rm(list=ls())

require(raster)
require(ncdf4)
require(xlsx)
require(gdalUtils)
###################################
###################################
## load all of the raw data sets ##
###################################
###################################

setwd("C:/Users/mloranty/Google Drive/GIS_Data/")

###################################################
### load GlobSnow SWE files for March and April ###
###################################################

##  ** NOTE TO HK FROM ML, 5/17 **##
# THIS IS MONTHLY DATA, INSTEAD OF DAILY/WEEKLY
# THAT WE ARE INTERESTED IN FOR SWE DEPLETION, 
# BUT THIS SCRIPT HAS GOOD EXAMPLES OF HOT TO LOAD 
# DESIRED NETCDF FILES AS RASTER OBJECTS
#
# HAVE A LOOK AND LET ME KNOW IF ANYTHING IS UNCLEAR

#list the monthly files
swe.files <- list.files(path='GlobSNOW/SWE/monthly_nc/',full.names=T)

# read data for March, 1980-2014
# make a Raster stack from a list of file paths - subset the filenames to include only months of interest
swe.march <- stack(swe.files[which(substr(swe.files,54,55)=='03')],varname='SWE_avg')
# EASE Grid Projection - this is not properly defined in the GlobSnow files #
# note there are some wonky projection issues with some of these data sets 
projection(swe.march) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"

# read data for April, 1980-2014
swe.april <- stack(swe.files[which(substr(swe.files,54,55)=='04')],varname='SWE_avg')
# EASE Grid Projection - this is not properly defined in the GlobSnow files #
projection(swe.april) <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"


## reclassify to make non-positive values NA
r <- matrix(c(-Inf,0,NA),ncol=3)

swe.april <- reclassify(swe.april,r,overwrite=T,
                        filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GlobSnow_SWE_April.tif')
swe.march <- reclassify(swe.march,r,overwrite=T,
                        filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GlobSnow_SWE_March.tif')

# calculate the difference FOR EACH YEAR 1980-2014
# this works on each layer in the stack since they have the same number of layers
swe.dif <- swe.march-swe.april
# reclassify to exclude positive, and negative data #
swe.dif.pos <- reclassify(swe.dif,matrix(c(-Inf,0,NA),ncol=3),
                          filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GlobSnow_SWE_Mar-Apr_pos.tif')
swe.dif.neg <- reclassify(swe.dif,matrix(c(0,Inf,NA),ncol=3),
                          filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GlobSnow_SWE_Mar-Apr_neg.tif')

## the following is not necessary for the swe decline, but has more examples


################################################
###  load CRU temp data for March and April  ###
################################################

cru.tmp <- brick('cru/cru_ts3.24.1901.2015.tmp.dat.nc')
cru.tmp <- crop(cru.tmp,c(-180,180,50,90))

# reproject to match GlobSnow #
cru.ease <- projectRaster(cru.tmp,swe.dif.pos,filename='cru/cru_ts3.24.1901.2015.tmp.dat.ease.N.tif',overwrite=T)

# subset March data for 1980-2014 #
cru.mar <- subset(cru.ease,
                  which(as.numeric(substr(cru.tmp@z$time,1,4))>1979 &
                   as.numeric(substr(cru.tmp@z$time,1,4))<2015 &
                   as.numeric(substr(cru.tmp@z$time,6,7))==3))

# subset April data for 1980-2014 #
cru.apr <- subset(cru.ease,
                  which(as.numeric(substr(cru.tmp@z$time,1,4))>1979 &
                          as.numeric(substr(cru.tmp@z$time,1,4))<2015 &
                          as.numeric(substr(cru.tmp@z$time,6,7))==4))

# subset May data for 1980-2014 #
cru.may <- subset(cru.ease,
                  which(as.numeric(substr(cru.tmp@z$time,1,4))>1979 &
                          as.numeric(substr(cru.tmp@z$time,1,4))<2015 &
                          as.numeric(substr(cru.tmp@z$time,6,7))==5))

writeRaster(cru.mar,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Mar.ease.N.tif')
writeRaster(cru.apr,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Apr.ease.N.tif')
writeRaster(cru.may,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.May.ease.N.tif')

# calculate the difference #
cru.dif <- cru.apr-cru.mar
writeRaster(cru.dif,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Mar_Apr.dif.ease.N.tif')

# reclass to include only positive values
cru.dif.pos <- reclassify(cru.dif,r,overwrite=T,
                          filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/cru_ts3.24.1980.2014.tmp.Apr_Mar_pos.ease.N.tif')

# calculate swe reduction per unit temp increase using positive datasets
swe.cru <- swe.dif.pos/cru.dif.pos
writeRaster(swe.cru,overwrite=T,
            filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/swe_cru_Mar_Apr_change.tif')

###################################################
### load, reproject, and mosaic MODIS VCF data  ###
###################################################
vcf.2014.w <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_west_mosaic_sin.tif')
vcf.2014.e <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_east_mosaic_sin.tif')
 
r <- matrix(c(101,255,NA),ncol=3)

vcf.2014.w <- reclassify(vcf.2014.w,r,progress='text', 
                          filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_west_mosaic_sin_rcl.tif')

vcf.2014.w.ease <- projectRaster(vcf.2014.w,cru.dif,progress='text',
                                 filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_west_mosaic_ease.tif')

vcf.2014.e <- reclassify(vcf.2014.e,r,progress='text',overwrite=T,  
                          filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_east_mosaic_sin_rcl.tif')

vcf.2014.e.ease <- projectRaster(vcf.2014.e,cru.dif,progress='text',overwrite=T,  
                                 filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_east_mosaic_ease.tif')

vcf.2014.ease <- mosaic(vcf.2014.w.ease,vcf.2014.e.ease,fun=mean,overwrite=T,  
                        filename='MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')

#vcf.2014.ease <- raster('MODIS/MOD44B/hdf_tiles/MOD44B_2014_mosaic_ease.tif')

###################################################
### load daily SWE files from Lawrence Mudryk    ##
### these are a combo of GlobSnow, MERRA & Brown ##
###################################################

###SKIP FOR NOW - HAVING PROBLEMS WITH EASE 2.0 GRID PROJECTION###
#setwd("Y:/swe_mudryk/")
swe <- brick('SWE_obsMEAN_M2BGS_ease2.2010.nc')

## note this is the correct PROJ.4 description for EASE2, but not sure this works with Lawrence's files ##
projection(swe) <- "+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m "
swe <- brick(swe)

sce <- stack('GlobSNOW/SCE/GlobSnow_SE_FSC_L3B-M_NH_201003_v2.1_1.nc',
             'GlobSNOW/SCE/GlobSnow_SE_FSC_L3B-M_NH_201004_v2.1_1.nc',
             'GlobSNOW/SCE/GlobSnow_SE_FSC_L3B-M_NH_201005_v2.1_1.nc',
             'GlobSNOW/SCE/GlobSnow_SE_FSC_L3B-M_NH_201006_v2.1_1.nc',varname='FSC')
  
# agb <- raster('boreal_biomass_thurner/thurner_biomass_v3_agb.tif') #thurner et al biomass
# agb.na <- raster('boreal_biomass_neigh/NACP_BOREAL_BIOME_BIOMASS_1273/data/NA_500_ecoLc_rcAGB.tif') # NASA biomass
# agb.eu <- raster('boreal_biomass_neigh/EURASIA_BIOME_1278/data/EU_500_ecoLc_rcagb.tif')

###################################################
### load GLC2000 landcover map                   ##
###                                              ##
###################################################
glc <- raster('GLC2000/Tiff/glc2000_v1_1.tif')
projection(glc) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
glc <- crop(glc,c(-180,180,50,90))
glc.legend <- read.xlsx('/Users/mloranty/Google Drive/GIS_Data/GLC2000/Tiff/Global_Legend.xls',sheetIndex=1)[,1:2]
colnames(glc.legend) <- c('zone','names')

#pr <- projectRaster(vcf.2014.ease,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",method='ngb')

## write function to calculate mode, for resampling
get.mode <- function(x,na.rm=T){
            f <- table(x)
            as.numeric(names(f)[which.max(f)])
}

## write function to calculate frequency of mode, for resampling
get.mode.freq <- function(x,na.rm=T){
                  f <- table(x)
                  (f)[which.max(f)]
}

## aggregate GLC2000 to ~0.25 degree resolution using mode, and make a map of the frequency
glc.mode <- aggregate(glc,fact=25,fun=get.mode)  
glc.mode.freq <- aggregate(glc,fact=25,fun=get.mode.freq) 

## reproject aggregated GLC2000 to EASE grid resolution
glc.ease <- projectRaster(glc.mode,cru.ease,method='ngb',progress='TEXT',overwrite=T,
                          filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode.tif')

glc.ease.freq <- projectRaster(glc.mode.freq,cru.ease,method='ngb',progress='TEXT',overwrite=T,
                          filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode_freq.tif')

## create mask files for 50% and 75% vegetation coverage in each class
glc.ease.50.mask <- reclassify(glc.ease.freq,matrix(c(0,312,NA,313,625,50),ncol=3,byrow=T),right=NA,overwrite=T,
                               filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode_freq_50.tif')

glc.ease.75.mask <- reclassify(glc.ease.freq,matrix(c(0,468,NA,469,615,75),ncol=3,byrow=T),right=NA,overwrite=T,
                               filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_mode_freq_75.tif')

## mask glc to include majority pixels only
glc.50pct.ease <- mask(glc.ease,glc.ease.50.mask,maskvalue=50,updatevalue=NA,inverse=T,overwrite=T,
                        filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_50_mask.tif')

glc.75pct.ease <- mask(glc.ease,glc.ease.75.mask,maskvalue=75,updatevalue=NA,inverse=T,overwrite=T,
                       filename='C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GLC2000_EASE_from_0.25_50_mask.tif')

## STOPPING HERE FOR NOW ##
dat <- ls()
save.image(dat,file="C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/monthly_swe_analysis_23Nov.RData")

#############################################
#############################################
##                                         ##
## LOAD ALL PREVIOUSLY PROCESSED DATA SETS ##
##                                         ##
#############################################
#############################################

swe.dif.pos <- raster('C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GlobSnow_SWE_Mar-Apr_pos.tif')
swe.dif.neg <- raster('C:/Users/mloranty/Google Drive/GIS_projects/boreal_snow_picker/data/GlobSnow_SWE_Mar-Apr_neg.tif')
#############################################
##      zONAL STATISTICS                   ##
#############################################

#make a table of cell counts 
t <- table(getValues(glc.ease))
t <- cbind(as.numeric(rownames(t)),as.numeric(t))
colnames(t) <- c('zone','cell.count')

#combine with glc legend
glc.table <- merge(glc.legend,t,all=T)
rm(t)

#calculate zonal statistics for each year
veg.stats <- zonal(swe.cru,glc.ease,fun='mean',na.rm=T)
colnames(veg.stats)[2:36] <- 1980:2014
rowMeans(veg.stats[,2:36],na.rm=T)

stf <- merge(glc.legend,zonal(swe.cru,glc.ease,fun='mean',na.rm=T))

stf <- cbind(stf,zonal(swe.tmp1,glc.laea,fun='sd',na.rm=T)[,2],
             zonal(swe.tmp2,glc.laea,fun='mean',na.rm=T)[,2],
             zonal(swe.tmp2,glc.laea,fun='sd',na.rm=T)[,2],
             zonal(swe.tmp3,glc.laea,fun='mean',na.rm=T)[,2],
             zonal(swe.tmp3,glc.laea,fun='sd',na.rm=T)[,2])

colnames(stf) <- c('zone','name','mean.MA','sd.MA','mean.AM','sd.AM','mean.MAM','sd.MAM')
stf[,3:8] <- round(stf[,3:8],digits=2)
write.xlsx(stf,file='GLC2000_stf_2010_zonal.xlsx')


agb.dnf <- mask(agb.laea,glc.laea,maskvalue=5,updatevalue=NA,inverse=T)
agb.enf <- calc(stack(mask(agb.laea,glc.laea,maskvalue=4,updatevalue=NA,inverse=T),
                      mask(agb.laea,glc.laea,maskvalue=9,updatevalue=NA,inverse=T)),fun=max)


vcf.enf <- calc(stack(mask(vcf.laea,glc.laea,maskvalue=4,updatevalue=NA,inverse=T),
                mask(vcf.laea,glc.laea,maskvalue=9,updatevalue=NA,inverse=T)),fun=max)
vcf.dnf <- mask(vcf.laea,glc.laea,maskvalue=5,updatevalue=NA,inverse=T)
vcf.mix <- mask(vcf.laea,glc.laea,maskvalue=6,updatevalue=NA,inverse=T)
vcf.dbf <- mask(vcf.laea,glc.laea,maskvalue=2,updatevalue=NA,inverse=T)

## preliminary plots
par(mfrow=c(2,2))
plot(vcf.2014.ease,swe.cru[[1]],ylim=c(-5,25),xlab='Tree Cover (%)',ylab='SWE/Temp',main="1980")
plot(vcf.2014.ease,swe.cru[[10]],ylim=c(-5,25),xlab='Tree Cover (%)',ylab='SWE/Temp',main="1990",col='red')
plot(vcf.2014.ease,swe.cru[[20]],ylim=c(-5,25),xlab='Tree Cover (%)',ylab='SWE/Temp',main="2000",col='blue')
plot(vcf.2014.ease,swe.cru[[30]],ylim=c(-5,25),xlab='Tree Cover (%)',ylab='SWE/Temp',main="2010",col='orange')


#############################################
#############################################
# OLD JUNK CODE
#############################################
#############################################
r.dnf <- lm(getValues(swe.tmp1)~getValues(vcf.dnf))
r.enf <- lm(getValues(swe.tmp1)~getValues(vcf.enf))
r.vcf <- lm(getValues(swe.tmp1)~getValues(vcf.laea))

plot(vcf.enf,swe.tmp1,ylim=c(-5,20))

plot(vcf.dnf,swe.tmp1,ylim=c(-5,10))


plot(vcf.laea,swe.tmp1,ylim=c(-5,15))

## make a plot using biomass data ##
pdf(6,10,file='march_april_agb_plots.pdf')
par(mfcol=c(3,1))
plot(agb.laea,swe.tmp1,pch=16,main='All Vegetation Types',cex=0.7,
     xlab=expression(paste('Biomass (kg C ',m^-2,')')),
     ylab=expression(paste('SWE Reduction Rate (mm ',C*degree^-1,')')),
     ylim=c(0,15),xlim=c(0,10))
r.agb <- lm(getValues(swe.tmp1)~getValues(agb.laea))
abline(r.agb,col='red',lwd=2)
legend('topright',paste('r2 =',round(summary(r.agb)$adj.r.squared,digits=3)),bty='n')

plot(agb.enf,swe.tmp1,pch=16,main='Evergreen Needleleaf',cex=0.7,
     xlab=expression(paste('Biomass (kg C ',m^-2,')')),
     ylab=expression(paste('SWE Reduction Rate (mm ',C*degree^-1,')')),
     ylim=c(0,15),xlim=c(0,10))
r.agb <- lm(getValues(swe.tmp1)~getValues(agb.enf))
abline(r.agb,col='red',lwd=2)
legend('topright',paste('r2 =',round(summary(r.agb)$adj.r.squared,digits=3)),bty='n')

plot(agb.dnf,swe.tmp1,pch=16,main='Deciduous Needleleaf',cex=0.7,
     xlab=expression(paste('Biomass (kg C ',m^-2,')')),
     ylab=expression(paste('SWE Reduction Rate (mm ',C*degree^-1,')')),
     ylim=c(0,15),xlim=c(0,10))
r.agb <- lm(getValues(swe.tmp1)~getValues(agb.dnf))
abline(r.agb,col='red',lwd=2)
legend('topright',paste('r2 =',round(summary(r.agb)$adj.r.squared,digits=3)),bty='n')
dev.off()
##
## make a plot using vcf data ##
pdf(6,10,file='march_april_vcf_plots.pdf')
par(mfcol=c(3,1))
plot(vcf.laea,swe.tmp1,pch=16,main='All Vegetation Types',cex=0.7,
     xlab='Canopy Cover (%)',
     ylab=expression(paste('SWE Reduction Rate (mm ',C*degree^-1,')')),
     ylim=c(0,15),xlim=c(0,100))
r.agb <- lm(getValues(swe.tmp1)~getValues(vcf.laea))
abline(r.agb,col='red',lwd=2)
legend('topleft',paste('r2 =',round(summary(r.agb)$adj.r.squared,digits=3)),bty='n')

plot(vcf.enf,swe.tmp1,pch=16,main='Evergreen Needleleaf',cex=0.7,
     xlab='Canopy Cover (%)',
     ylab=expression(paste('SWE Reduction Rate (mm ',C*degree^-1,')')),
     ylim=c(0,15),xlim=c(0,100))
r.agb <- lm(getValues(swe.tmp1)~getValues(vcf.enf))
abline(r.agb,col='red',lwd=2)
legend('topleft',paste('r2 =',round(summary(r.agb)$adj.r.squared,digits=3)),bty='n')

plot(vcf.dnf,swe.tmp1,pch=16,main='Deciduous Needleleaf',cex=0.7,
     xlab='Canopy Cover (%)',
     ylab=expression(paste('SWE Reduction Rate (mm ',C*degree^-1,')')),
     ylim=c(0,15),xlim=c(0,100))
r.agb <- lm(getValues(swe.tmp1)~getValues(vcf.dnf))
abline(r.agb,col='red',lwd=2)
legend('topleft',paste('r2 =',round(summary(r.agb)$adj.r.squared,digits=3)),bty='n')
dev.off()
##