#############
# LIBRARIES #
#############
.libPaths(c('/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library/',.libPaths()))
library(celestial)
library(devtools)
#
library(Cairo)
library(ProFound) 
library(magicaxis)
library(data.table)
library(Rcpp)
require(foreign)
require(MASS)
library(Rfits)
#
###################
# Setup PATHNAMES #
###################
#
inputargs=commandArgs(TRUE)
ra=as.character(inputargs[1])
dec=as.character(inputargs[2])

InputDir="/group/pawsey0160/waves/profound/tiles_0.3Res/"
OutputDir="/group/pawsey0160/waves/profound/detect_0.3Res/"
PlotDir="/group/pawsey0160/waves/profound/plots/"
#
#################################
# reading in the SWarped images #
#################################

print(paste0('Reading in SWarped images for: ', ra, '_', dec))

kids_r=Rfits_read_image(paste0(InputDir, 'r_', ra, '_', dec, '.fits'), ext=1)
kids_i=Rfits_read_image(paste0(InputDir, 'i_', ra, '_', dec, '.fits'), ext=1)
viking_Z=Rfits_read_image(paste0(InputDir, 'Z_', ra, '_', dec, '.fits'), ext=1)
viking_Y=Rfits_read_image(paste0(InputDir, 'Y_', ra, '_', dec, '.fits'), ext=1)

############################
# generate detection image #
############################

InputImages=list(kids_r, kids_i, viking_Z, viking_Y)
pixcut_value=9
skycut_value=1.5

print(paste0('Starting detection stack for: ', ra, '_', dec))


# running profound with a high skyvut in order to define mask. 
# Skycut value is important here. A higher value will mean that fewer noisy regions are registered as objects, which is good, as they will
# then be picked up as a high sky, which can then be masked. 
pro_r_test = profoundProFound(kids_r, skycut=4, roughpedestal=T, grid=50, boxiters=4, verbose=F, SBdilate=2, pixcut = pixcut_value, tolerance = 15, reltol = -10)
pro_i_test = profoundProFound(kids_i, skycut=4, roughpedestal=T, grid=50, boxiters=4, verbose=F, SBdilate=2, pixcut = pixcut_value, tolerance = 15, reltol = -10)
pro_Z_test = profoundProFound(viking_Z, skycut=10, roughpedestal=T, grid=50, boxiters=4, verbose=F, SBdilate=2, pixcut = pixcut_value, tolerance = 15, reltol = -10)
pro_Y_test = profoundProFound(viking_Y, skycut=10, roughpedestal=T, grid=50, boxiters=4, verbose=F, SBdilate=2, pixcut = pixcut_value, tolerance = 15, reltol = -10)

print(paste0('Making initial sky grid for: ', ra, '_', dec))
sky_r_test = profoundMakeSkyGrid(kids_r$imDat, objects=pro_r_test$objects_redo, doclip=T, grid=50, box=100)
sky_i_test = profoundMakeSkyGrid(kids_i$imDat, objects=pro_i_test$objects_redo, doclip=T, grid=50, box=100)
sky_Z_test = profoundMakeSkyGrid(viking_Z$imDat, objects=pro_Z_test$objects_redo, doclip=T, grid=50, box=100)
sky_Y_test = profoundMakeSkyGrid(viking_Y$imDat, objects=pro_Y_test$objects_redo, doclip=T, grid=50, box=100)

labelCoord=c(40, 950)
png(filename=paste0(PlotDir,'SkyRMS_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
par(mfrow=c(2, 2),mar=c(2,2,2,2))
magimage(sky_r_test$skyRMS, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
magimage(sky_i_test$skyRMS, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
magimage(sky_Z_test$skyRMS, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
magimage(sky_Y_test$skyRMS, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
dev.off()

labelCoord=c(40, 950)
png(filename=paste0(PlotDir,'Sky_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
par(mfrow=c(2, 2),mar=c(2,2,2,2))
magimage(sky_r_test$sky, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
magimage(sky_i_test$sky, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
magimage(sky_Z_test$sky, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
magimage(sky_Y_test$sky, bad=0.5)
text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
dev.off()


# based on the above sky measurement in each band, here we select which pixels in each detection band to mask. The higher the modifier value, the more
# pixels become masked. 
# also making sure that the zero pixels are masked
mask_r_test = (sky_r_test$sky > (2*median(sky_r_test$sky,na.rm = T) - 0.8*quantile(sky_r_test$sky,0.01,na.rm=TRUE)) )
mask_i_test = (sky_i_test$sky > (2*median(sky_i_test$sky,na.rm = T) - 2.0*quantile(sky_i_test$sky,0.01,na.rm=TRUE)) )
mask_Z_test = (sky_Z_test$skyRMS > (2*median(sky_Z_test$skyRMS,na.rm = T) - 1.0*quantile(sky_Z_test$skyRMS,0.01,na.rm=TRUE)) )
mask_Y_test = (sky_Y_test$skyRMS > (2*median(sky_Y_test$skyRMS,na.rm = T) - 1.0*quantile(sky_Y_test$skyRMS,0.01,na.rm=TRUE)) )


# this makes sure that we're not masking actual astrophysical objects. 
# The higher the modifier value, 
BlurIm=profoundImBlur(kids_r$imDat,25)
unmask_test = BlurIm > 1.2e-11
mask_r_test[unmask_test] = 0L
mask_i_test[unmask_test] = 0L
mask_Z_test[unmask_test] = 0L
mask_Y_test[unmask_test] = 0L

png(filename=paste0(PlotDir,'Unmask_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
magimage(BlurIm); magimage(unmask_test, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
dev.off()

# now to mask all regions that don't have data:
mask_r_test[kids_r$imDat == 0] = 1L
mask_i_test[kids_i$imDat == 0] = 1L
mask_Z_test[viking_Z$imDat == 0] = 1L
mask_Y_test[viking_Y$imDat == 0] = 1L


# now checking for any pixels that are masked in all four bands
AllMask = (mask_r_test==1) & (mask_i_test==1) & (mask_Z_test==1) & (mask_Y_test==1) & (kids_i$imDat!=0) # only specifying pixels with i-band imaging. 
mask_i_test[AllMask] = 0L # ensuring they're at least unmasked in the i-band, which seems the least susceptible to noise features. 


labelCoord=c(40, 950)
png(filename=paste0(PlotDir,'DetectBands_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
par(mfrow=c(2, 2),mar=c(2,2,2,2))
magimage(pro_r_test$image, bad=0.5); magimage(mask_r_test, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
magimage(pro_i_test$image, bad=0.5); magimage(mask_i_test, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
magimage(pro_Z_test$image, bad=0.5); magimage(mask_Z_test, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
magimage(pro_Y_test$image, bad=0.5); magimage(mask_Y_test, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
dev.off()

stack_test = profoundMakeStack(image_list = list(pro_r_test$image, pro_i_test$image, pro_Z_test$image, pro_Y_test$image), 
                         sky_list = list(pro_r_test$sky, pro_i_test$sky, pro_Z_test$sky, pro_Y_test$sky), 
                         skyRMS_list = list(pro_r_test$skyRMS, pro_i_test$skyRMS, pro_Z_test$skyRMS, pro_Y_test$skyRMS),
                         mask_list = list(mask_r_test, mask_i_test, mask_Z_test, mask_Y_test),
                         magzero_in = c(0,0,30,30)
                         )

stack_sky = profoundMakeStack(image_list = list(pro_r_test$sky, pro_i_test$sky, pro_Z_test$sky, pro_Y_test$sky), 
                         skyRMS_list = list(pro_r_test$skyRMS, pro_i_test$skyRMS, pro_Z_test$skyRMS, pro_Y_test$skyRMS),
                         mask_list = list(mask_r_test, mask_i_test, mask_Z_test, mask_Y_test),
                         magzero_in = c(0,0,30,30)
                         )




png(filename=paste0(PlotDir,'DetectStack_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
magimage(stack_test$image, bad=0.5)
dev.off()

print(paste0('Running proFound in detect mode for: ', ra, '_', dec))

#######################
# Run ProFound Script #
#######################
#
# needs to be hdr not header. 
everything=profoundProFound(image=stack_test$image+stack_sky$image, 
							header=kids_r$hdr, 
							sky=stack_sky$image, 
							skycut = skycut_value, 
							pixcut = pixcut_value, 
							ext = 1, 
							tolerance = 15, 
							reltol = -10, 
							keepsegims = TRUE, 
							magzero = stack_test$magzero, 
							dotot = TRUE, 
							docol = TRUE, 
							dogrp = FALSE, 
							verbose = TRUE, 
							boxiters = 4, 
							box=100, 
							grid = 50, 
							roughpedestal = TRUE, 
							stats = TRUE, 
							groupstats = TRUE, 
							mask = 0, 
							appdiam = 1.4, 
							fluxtype = 'Jansky', 
							SBdilate=2, 
							redosky=FALSE) 

print(paste0('Saving outputs for: ', ra, '_', dec))

# adding the detection masks into the output
everything$detect_mask_r=mask_r_test
everything$detect_mask_i=mask_i_test
everything$detect_mask_Z=mask_Z_test
everything$detect_mask_Y=mask_Y_test


cat('---saving segim plot---\n')
png(file = paste0(PlotDir,'segimPlot_', ra, '_', dec, '.png'), width = 1000, height = 1000)
profoundSegimPlot(everything)
dev.off()

# cat('---saving RDS---\n')
# saveRDS(everything,file=paste0(OutputDir,'waves_detect_', ra, '_', dec, '.rds'))

cat('---saving diagnostics plot---\n')
CairoPDF(file=paste0(PlotDir,'diagnostics_', ra, '_', dec, '.pdf'),width=24.0,height=24.0)
plot(everything)
dev.off()

# making sure that the groupID values and Ngroup values are in the groupstats and segstats data frames. 
temp=data.frame(segID=unlist(everything$group$groupsegID$segID), 
					groupID=rep(everything$group$groupsegID$groupID, times=everything$group$groupsegID$Ngroup))

everything$groupstats$Ngroup=everything$group$groupsegID$Ngroup[match(everything$groupstats$groupID, everything$group$groupsegID$groupID)]
everything$segstats$groupID=temp$groupID[match(everything$segstats$segID, temp$segID)]

cat('---saving FITS---\n')
Rfits_write(everything, filename=paste0(OutputDir,'waves_detect_', ra, '_', dec, '.fits'), compress=TRUE, list_sub=c('segim', 'segim_orig', 'objects', 'objects_redo', 'sky', 
			'skyRMS', 'image', 'mask', 'segstats', 'groupstats', 'detect_mask_r', 'detect_mask_i', 'detect_mask_Z', 'detect_mask_Y'))

Rfits_write_header(filename=paste0(OutputDir,'waves_detect_', ra, '_', dec, '.fits'), ext=1, 
				   keyvalues=list(Nseg=everything$Nseg, 
				   magzero=everything$magzero, 
				   pixscale=everything$pixscale, 
				   imarea=everything$imarea, 
				   skyLL=everything$skyLL, 
				   skyChiSq=everything$skyChiSq, 
				   time=everything$time,
				   R_version=everything$R.version$version.string, 
				   ProFound_version=everything$ProFound.version))
