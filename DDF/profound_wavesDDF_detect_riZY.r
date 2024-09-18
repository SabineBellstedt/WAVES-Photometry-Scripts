# remaking the detection script for the DDF fields using a riZY stack

library(celestial)
library(devtools)
#
library(Cairo)
library(ProFound) 
library(ProPane)
library(magicaxis)
library(data.table)
library(Rcpp)
require(foreign)
require(MASS)
library(Rfits)
library(foreach)
library(doParallel)
CoreNumber = 7


InputDir="/Volumes/WAVES/waves/profound/ddfProMosv2/"
InputDir2="/Volumes/WAVES/waves/profound/ddfProMosv3/"
OutputDir = '/Volumes/ThunderBay/DDF/profound_ddf_v1/detect/'
PlotDir="/Volumes/ThunderBay/DDF/profound_ddf_v1/plots/"

TileCat=fread('/Volumes/ThunderBay/DDF/profound_ddf_v1/wd10tiles_names.csv')


##################################
# reading in the ProPaned images #
##################################

# looping through each of the sq degrees
registerDoParallel(cores=CoreNumber)

# foreach(ii=(1:length(TileCat$ra))[-c(2)])%dopar%{
foreach(ii=(1:length(TileCat$ra)))%dopar%{
# for(ii in 1:length(TileCat$ra)){
# for(ii in 2:2){
# for(ii in 55:60){
	ra = format(round(TileCat$ra[ii], 1), nsmall = 1)
    dec = format(round(TileCat$dec[ii], 1), nsmall = 1)

	print(paste0('Starting process for ', ra, '_', dec))
	
	StartTime = Sys.time()

	# only do this process of the detect output doesn't already exist
	if(!file.exists(paste0(OutputDir,'waves_detect_', ra, '_', dec, '.fits'))){
		zz <- file(paste0(OutputDir,'waves_detect_log_', ra, '_', dec, '.txt'), open = "wt")
		sink(zz, type='message')


		message(paste0('Reading in ProMo images for: ', ra, '_', dec))
	
		Available_Images=c(FALSE, FALSE, FALSE, FALSE)
		
		if(file.exists(paste0(InputDir, 'r_', ra, '_', dec, '.fits'))){
			vst_r=Rfits_read(paste0(InputDir, 'r_', ra, '_', dec, '.fits'), pointer=FALSE) # includes $image $weight and $inVar
			vst_r$image$imDat[!is.finite(vst_r$image$imDat)] = 0L
			Available_Images[1]=TRUE
		}
		if(file.exists(paste0(InputDir, 'i_', ra, '_', dec, '.fits'))){
			vst_i=Rfits_read(paste0(InputDir, 'i_', ra, '_', dec, '.fits'), pointer=FALSE)
			vst_i$image$imDat[!is.finite(vst_i$image$imDat)] = 0L
			Available_Images[2]=TRUE
		}
		if(file.exists(paste0(InputDir2, 'Z_', ra, '_', dec, '.fits'))){
			viking_Z=Rfits_read(paste0(InputDir2, 'Z_', ra, '_', dec, '.fits'), pointer=FALSE) 
			viking_Z$image$imDat[!is.finite(viking_Z$image$imDat)] = 0L
			Available_Images[3]=TRUE
		}
		if(file.exists(paste0(InputDir2, 'Y_', ra, '_', dec, '.fits'))){
			viking_Y=Rfits_read(paste0(InputDir2, 'Y_', ra, '_', dec, '.fits'), pointer=FALSE)
			viking_Y$image$imDat[!is.finite(viking_Y$image$imDat)] = 0L
			Available_Images[4]=TRUE
		}
		
		if(sum(Available_Images)>0){ # only continue if we have any images available
		
			# for those images that are missing, replicate the image in an existing band
			index=1:4
			FirstAvailableImage=index[Available_Images>0][1]
			if(FirstAvailableImage==1){
				DummyImage = vst_r
			}else if(FirstAvailableImage==2){
				DummyImage = vst_i
			}else if(FirstAvailableImage==3){
				DummyImage = viking_Z
			}else if(FirstAvailableImage==4){
				DummyImage = viking_Y
			}
		
			if(Available_Images[1]==FALSE){
				vst_r = DummyImage
				vst_r$image$imDat[]=NA
				vst_r$weight$imDat[]=NA
				vst_r$inVar$imDat[]=NA
				message('*** Missing image in the r band. Generating all-NA image ***')
			}
			if(Available_Images[2]==FALSE){
				vst_i = DummyImage
				vst_i$image$imDat[]=NA
				vst_i$weight$imDat[]=NA
				vst_i$inVar$imDat[]=NA
				message('*** Missing image in the i band. Generating all-NA image ***')
			}
			if(Available_Images[3]==FALSE){
				viking_Z = DummyImage
				viking_Z$image$imDat[]=NA
				viking_Z$weight$imDat[]=NA
				viking_Z$inVar$imDat[]=NA
				message('*** Missing image in the Z band. Generating all-NA image ***')
			}
			if(Available_Images[4]==FALSE){
				viking_Y = DummyImage
				viking_Y$image$imDat[]=NA
				viking_Y$weight$imDat[]=NA
				viking_Y$inVar$imDat[]=NA
				message('*** Missing image in the Y band. Generating all-NA image ***')
			}

			labelCoord=c(40, 950)
			CairoPNG(filename=paste0(PlotDir,'Tiles_qdiff_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			par(mfrow=c(2, 2),mar=c(2,2,2,2))
			magimage(vst_r$image$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
			magimage(vst_i$image$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
			magimage(viking_Z$image$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
			magimage(viking_Y$image$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
			dev.off()

			CairoPNG(filename=paste0(PlotDir,'Tiles_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			par(mfrow=c(2, 2),mar=c(2,2,2,2))
			magimage(vst_r$image$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
			magimage(vst_i$image$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
			magimage(viking_Z$image$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
			magimage(viking_Y$image$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
			dev.off()

			# actually need to mask out all pixels in the image with weight==0
			vst_r$image$imDat[vst_r$weight$imDat==0] = NA
			vst_i$image$imDat[vst_i$weight$imDat==0] = NA
			viking_Z$image$imDat[viking_Z$weight$imDat==0] = NA
			viking_Y$image$imDat[viking_Y$weight$imDat==0] = NA		

			Weight_scales = list(r=median(vst_r$inVar$imDat/vst_r$weight$imDat, na.rm=TRUE), 
				i= median(vst_i$inVar$imDat/vst_i$weight$imDat, na.rm=TRUE), 
				Z=median(viking_Z$inVar$imDat/viking_Z$weight$imDat, na.rm=TRUE), 
				Y=median(viking_Y$inVar$imDat/viking_Y$weight$imDat, na.rm=TRUE))	

			# run imaging through propaneBadPix to remove any bad pixels. 
			message(paste0('Conducting ProPane bad pixel removal for : ', ra, '_', dec))

			vst_r$image_badPix = propaneBadPix(vst_r$image, inVar=vst_r$weight$imDat*Weight_scales$r, sigma = 20, pixcut=1, dilate=TRUE, size=5)
			vst_i$image_badPix = propaneBadPix(vst_i$image, inVar=vst_i$weight$imDat*Weight_scales$i, sigma = 20, pixcut=1, dilate=TRUE, size=5)
			viking_Z$image_badPix = propaneBadPix(viking_Z$image, inVar=viking_Z$weight$imDat*Weight_scales$Z, sigma=20, pixcut=1, dilate=TRUE, size=5, cold=TRUE) # cold=TRUE removes very -ve pixels, which covers snowflaking
			viking_Y$image_badPix = propaneBadPix(viking_Y$image, inVar=viking_Y$weight$imDat*Weight_scales$Y, sigma=20, pixcut=1, dilate=TRUE, size=5, cold=TRUE)

			CairoPNG(filename=paste0(PlotDir,'Tiles_badPixRemoved_qdiff_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			par(mfrow=c(2, 2),mar=c(2,2,2,2))
			magimage(vst_r$image_badPix$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
			magimage(vst_i$image_badPix$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
			magimage(viking_Z$image_badPix$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
			magimage(viking_Y$image_badPix$imDat, bad=0.5,qdiff=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
			dev.off()

			CairoPNG(filename=paste0(PlotDir,'Tiles_badPixRemoved_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			par(mfrow=c(2, 2),mar=c(2,2,2,2))
			magimage(vst_r$image_badPix$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
			magimage(vst_i$image_badPix$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
			magimage(viking_Z$image_badPix$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
			magimage(viking_Y$image_badPix$imDat, bad=0.5)
			text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
			dev.off()

			pixcut_value = 9
			skycut_value = 1.5 # for detect run
			
			magzero_VST = 23.9 # image units are in microJy
			magzero_Viking = 23.9 # image units are in microJy
			
			message(paste0('Time taken since script beginning: ', Sys.time() - StartTime))
			message(paste0('Starting detection stack for: ', ra, '_', dec))
			
			

			


			# maybe we don't actually bother with the masking at all to make the stack, and then just use the inVar for the skyRMS into ProFOund

			# now to mask all regions that don't have data:
			mask_r = vst_r$weight$imDat == 0
			mask_i = vst_i$weight$imDat == 0
			mask_Z = viking_Z$weight$imDat == 0
			mask_Y = viking_Y$weight$imDat == 0

			# masking any obscenely hot pixels
			mask_r[vst_r$image$imDat > 1e10] == 1L
			mask_i[vst_i$image$imDat > 1e10] == 1L
			mask_Z[viking_Z$image$imDat > 1e10] == 1L
			mask_Y[viking_Y$image$imDat > 1e10] == 1L

			labelCoord=c(40, 950)
			CairoPNG(filename=paste0(PlotDir,'DetectBands_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			par(mfrow=c(2, 2),mar=c(2,2,2,2))
			magimage(vst_r$image_badPix$imDat, bad=0.5); magimage(mask_r, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
			magimage(vst_i$image_badPix$imDat, bad=0.5); magimage(mask_i, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
			magimage(viking_Z$image_badPix$imDat, bad=0.5); magimage(mask_Z, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
			magimage(viking_Y$image_badPix$imDat, bad=0.5); magimage(mask_Y, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
			dev.off()

			labelCoord=c(40, 950)
			CairoPNG(filename=paste0(PlotDir,'DetectBands_qdiff_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			par(mfrow=c(2, 2),mar=c(2,2,2,2))
			magimage(vst_r$image_badPix$imDat, bad=0.5, qdiff=TRUE, stretchscale=1960); magimage(mask_r, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
			magimage(vst_i$image_badPix$imDat, bad=0.5, qdiff=TRUE, stretchscale=1960); magimage(mask_i, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
			magimage(viking_Z$image_badPix$imDat, bad=0.5, qdiff=TRUE, stretchscale=1960); magimage(mask_Z, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
			magimage(viking_Y$image_badPix$imDat, bad=0.5, qdiff=TRUE, stretchscale=1960); magimage(mask_Y, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
			dev.off()


				

			skyRMS_list = list(
				1/sqrt(vst_r$weight$imDat*Weight_scales$r), 
				1/sqrt(vst_i$weight$imDat*Weight_scales$i), 
				1/sqrt(viking_Z$weight$imDat*Weight_scales$Z), 
				1/sqrt(viking_Y$weight$imDat*Weight_scales$Y))

			skyRMS_list[[1]][!is.finite(skyRMS_list[[1]])] = NA
			skyRMS_list[[2]][!is.finite(skyRMS_list[[2]])] = NA
			skyRMS_list[[3]][!is.finite(skyRMS_list[[3]])] = NA
			skyRMS_list[[4]][!is.finite(skyRMS_list[[4]])] = NA


			CairoPNG(filename=paste0(PlotDir,'skyRMS_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			par(mfrow=c(2, 2),mar=c(2,2,2,2))
			magimage(skyRMS_list[[1]], bad=0.5); magimage(mask_r, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="r",col="green",pos=4,cex=3)
			magimage(skyRMS_list[[2]], bad=0.5); magimage(mask_i, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="i",col="green",pos=4,cex=3)
			magimage(skyRMS_list[[3]], bad=0.5); magimage(mask_Z, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Z",col="green",pos=4,cex=3)
			magimage(skyRMS_list[[4]], bad=0.5); magimage(mask_Y, magmap=F, col=c(NA,hsv(alpha=0.3)),add=TRUE)
			text(x=labelCoord[1],y=labelCoord[2],label="Y",col="green",pos=4,cex=3)
			dev.off()

			# now patching any bad pixels
			message(paste0('Conducting ProPane pixel patching for : ', ra, '_', dec))
			vst_r$image_badPix = propanePatchPix(vst_r$image_badPix)
			vst_i$image_badPix = propanePatchPix(vst_i$image_badPix)
			viking_Z$image_badPix = propanePatchPix(viking_Z$image_badPix)
			viking_Y$image_badPix = propanePatchPix(viking_Y$image_badPix)


			message(paste0('Time taken since script beginning: ', Sys.time() - StartTime))
			message('Making stacked image')
			stack_detect = propaneStackFlatInVar(image_list = list(vst_r$image_badPix$imDat, vst_i$image_badPix$imDat, viking_Z$image_badPix$imDat, viking_Y$image_badPix$imDat), 
			                         # sky_list = list(pro_Z_disco$sky, pro_Y_disco$sky, pro_J_disco$sky, pro_H_disco$sky), 
			                         skyRMS_list = skyRMS_list,
			                         mask_list = list(mask_r, mask_i, mask_Z, mask_Y),
			                         magzero_in = c(magzero_VST,magzero_VST,magzero_Viking,magzero_Viking)
			                         )

			# stack_median = propaneStackFlatFunc(image_list = list(vst_r$image_badPix$imDat, vst_i$image_badPix$imDat, viking_Z$image_badPix$imDat, viking_Y$image_badPix$imDat))
			
			message(paste0('Time taken since script beginning: ', Sys.time() - StartTime))

			CairoPNG(filename=paste0(PlotDir,'DetectStack_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			magimage(stack_detect$image, bad=0.5)
			dev.off()

			CairoPNG(filename=paste0(PlotDir,'DetectStack_qdiff_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			magimage(stack_detect$image, bad=0.5, qdiff=TRUE)
			dev.off()

			CairoPNG(filename=paste0(PlotDir,'DetectStack_skyRMS_', ra, '_', dec, '.png'),width=40.0,height=40.0,units="cm",bg="transparent",res=240)
			magimage(stack_detect$skyRMS, bad=0.5)
			dev.off()

			message(paste0('Running proFound in detect mode for: ', ra, '_', dec))
			
			#######################
			# Run ProFound Script #
			#######################
			#
			# needs to be hdr not header. 
			everything=profoundProFound(image=stack_detect$image,#+stack_sky$image, 
										header=viking_Z$image$hdr, 
										# sky=stack_sky$image, 
										skyRMS=stack_detect$skyRMS,
										skycut = skycut_value, 
										pixcut = pixcut_value, 
										ext = 1, 
										tolerance = 15, 
										reltol = -10, 
										keepsegims = TRUE, 
										magzero = stack_detect$magzero, 
										dotot = TRUE, 
										docol = TRUE, 
										dogrp = FALSE, 
										verbose = TRUE, 
										boxiters = 4, 
										box=100, 
										grid = 50, 
										roughpedestal = FALSE, 
										stats = TRUE, 
										groupstats = TRUE, 
										mask = 0, 
										appdiam = 1.4, 
										fluxtype = 'Jansky', 
										SBdilate=2, 
										redosky=FALSE) 

			print(paste0('Saving outputs for: ', ra, '_', dec))
			
			# # adding the detection masH into the output
			everything$detect_mask_r=mask_r
			everything$detect_mask_i=mask_i
			everything$detect_mask_Z=mask_Z
			everything$detect_mask_Y=mask_Y

			everything$weight_scales=Weight_scales
			
			message('---saving segim plot---\n')
			CairoPNG(file = paste0(PlotDir,'segimPlot_', ra, '_', dec, '.png'), width = 1000, height = 1000)
			profoundSegimPlot(everything)
			dev.off()
			
			cat('---saving RDS---\n')
			saveRDS(everything,file=paste0(OutputDir,'waves_detect_', ra, '_', dec, '.rds'))
			
			
			message('---saving diagnostics plot---\n')
			# pdf(file=paste0(PlotDir,'diagnostics_', ra, '_', dec, '.pdf'),width=24.0,height=24.0)
			# plot(everything)
			# dev.off()
			CairoPNG(file=paste0(PlotDir,'diagnostics_', ra, '_', dec, '.png'),width = 2000, height = 2000)
			plot(everything)
			dev.off()
			
			# making sure that the groupID values and Ngroup values are in the groupstats and segstats data frames. 
			temp=data.frame(segID=unlist(everything$group$groupsegID$segID), 
								groupID=rep(everything$group$groupsegID$groupID, times=everything$group$groupsegID$Ngroup))
			
			everything$groupstats$Ngroup=everything$group$groupsegID$Ngroup[match(everything$groupstats$groupID, everything$group$groupsegID$groupID)]
			everything$segstats$groupID=temp$groupID[match(everything$segstats$segID, temp$segID)]
			
			message('---saving FITS---\n')
			Rfits_write(everything, filename=paste0(OutputDir,'waves_detect_', ra, '_', dec, '.fits'), compress=TRUE, 
				list_sub=c('segim', 'segim_orig', 'objects', 'objects_redo', 'sky', 
						'skyRMS', 'image', 'mask', 'segstats', 'groupstats', 'detect_mask_r', 'detect_mask_i', 'detect_mask_Z', 'detect_mask_Y', 'weight_scales'))
			
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

			sink()
			
		}else{
			message('no images available in this tile')
		}

	
	
	}else{
		print(paste0(OutputDir,'waves_detect_', ra, '_', dec, '.fits already exists'))
	}
}	



