library(celestial)
#
library(ProFound) 
library(magicaxis)
library(data.table)
library(Rcpp)
require(foreign)
require(MASS)
library(Rfits)
library(plotrix)
library(Rwcs)
library(doParallel)
library(bit64)
CoreNumber = 9

#
#  DEFINE FUNCTIONs
#
magAB2Jansky=function(x){10^(-0.4*(x-8.9))} # taken from ProSpect
#


autoDisco = function(skyRMS, blurSigma=200, kernelSmoothing=0.5, blur=TRUE){
	if(blur){
		blurSkyRMS = profoundImBlur(skyRMS,blurSigma)
		density_out = density(blurSkyRMS, bw=kernelSmoothing)
	}else{
		blurSkyRMS = skyRMS
		density_out = density(skyRMS, bw=kernelSmoothing)
	}
	
	Troughs=c()
	for(ii in 2:(length(density_out$x)-1)){
		if((density_out$y[ii+1] > density_out$y[ii]) & (density_out$y[ii-1] > density_out$y[ii])){
			Troughs=c(Troughs, density_out$x[ii])
		}
	}
	Troughs=c(Troughs, max(density_out$x))

	disco_reg = skyRMS*0
	for(value in 1:length(Troughs)){
		if(value==1){
			Mask = (blurSkyRMS<Troughs[value])
		}else{
			Mask = (blurSkyRMS<Troughs[value]) & (blurSkyRMS>=Troughs[value-1])
		}
		disco_reg[Mask] = value
	}	
	return(disco_reg)
}

################
#
# MAIN CODE
#
################
DetectDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/detect/'
MeasureDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/measure_WISE/'
PostprocessDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/postprocess_WISE_zpTest_GAAPoffset/'
PlotDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/plots/'
RefDir="/Volumes/ThunderBay/DDF/profound_ddf_v0/ref_Sabine/"

InputTargetCat=fread('/Volumes/ThunderBay/DDF/profound_ddf_v1/ddftiles_names.csv')
# InputTargetCat=fread('/Volumes/ThunderBay/DDF/profound_ddf_v1/ddftiles-2-3-10_names.csv')

# adding corrections for the ZP offsets based on colour distributions
ZP_u = list('WD01' = -0.1799418 + -0.07617032, 'WD02' = -0.08363246 + 0.02220393, 'WD03' = -0.1667924 + -0.06159079,  'WD10' = -0.07139491 + -0.003796734)
ZP_g = list('WD01' = -0.2216751 + -0.01626998, 'WD02' = -0.1296248  + 0.01620855,  'WD03' = -0.1237704 + 0.01272416,  'WD10' = -0.1243195 + -0.01465413)
ZP_r = list('WD01' = -0.2672564 + 0.002462328, 'WD02' = -0.09285457 + -0.003867115, 'WD03' = -0.06512765 + -0.02223372, 'WD10' = -0.06298653 + -0.02525463)
ZP_i = list('WD01' = -0.1177713 + -0.013364, 'WD02' =  0.01909506 + 0.0009301004,  'WD03' = -0.02955972 + -0.02821662, 'WD10' = 0.03467224 + -0.02472555)
ZP_Z = list('WD01' = -0.01062275 + 0.001333254, 'WD02' = -0.01024909 + -0.01010921,'WD03' = -0.02754233 + -0.01247592, 'WD10' = -0.01912775 + -0.02154232)

# compiling the above into one master list
ZP_corrections = list('u' = ZP_u, 
	'g' = ZP_g, 
	'r' = ZP_r, 
	'i' = ZP_i, 
	'Z' = ZP_Z)

# WD01
# u 0.0111292
# g 0.0763592
# r 0.158003
# i (+0.111071 + +0.109742)/2
# WD02
# u 0.0451242
# g -0.0403792
# r -0.0239738
# i (-0.0997323 + -0.103831)/2
# Z -0.0313122
# Y 0.00499546
# J -0.0593764
# H -0.0201537
# K -0.0172672
# WD03
# u 0.05033
# g 0.0178859
# r 0.0025961
# i (-0.0097085 + -0.00548141)/2
# Z 0.0270984
# Y -0.00827288
# J -0.0191033
# H 0.0100744
# K 0.00181977
# WD10
# u 0.0515447
# g -0.024937
# r -0.0158557
# i (-0.0617765 + -0.0663206)/2
# Z 0.0240375
# Y 0.041498
# J -0.00417184
# H 0.0113635
# K 0.0145904

# now instead computing ZP offsets based on the GAAP photometry
ZP_u = list('WD01' = -0.0111292, 'WD02' = -0.0451242, 'WD03' = -0.05033, 'WD10' = -0.0515447)
ZP_g = list('WD01' = -0.0763592, 'WD02' = +0.0403792, 'WD03' = -0.0178859, 'WD10' = +0.024937)
ZP_r = list('WD01' = -0.158003, 'WD02' = +0.0239738, 'WD03' = -0.0025961, 'WD10' = +0.0158557)
ZP_i = list('WD01' = (-0.111071 + -0.109742)/2, 'WD02' = (+0.0997323 + +0.103831)/2, 'WD03' = (+0.0097085 + +0.00548141)/2, 'WD10' = (+0.0617765 + +0.0663206)/2)
ZP_Z = list('WD01' = 0.0, 'WD02' = +0.0313122, 'WD03' = -0.0270984, 'WD10' = -0.0240375)
ZP_Y = list('WD01' = 0.0, 'WD02' = -0.00499546, 'WD03' = +0.00827288, 'WD10' = -0.041498)
ZP_J = list('WD01' = 0.0, 'WD02' = +0.0593764, 'WD03' = +0.0191033, 'WD10' = +0.00417184)
ZP_H = list('WD01' = 0.0, 'WD02' = +0.0201537, 'WD03' = -0.0100744, 'WD10' = -0.0113635)
ZP_K = list('WD01' = 0.0, 'WD02' = +0.0172672, 'WD03' = -0.00181977, 'WD10' = -0.0145904)



# compiling the above into one master list
ZP_corrections = list('u' = ZP_u, 
	'g' = ZP_g, 
	'r' = ZP_r, 
	'i' = ZP_i, 
	'Z' = ZP_Z, 
	'Y' = ZP_Y,
	'J' = ZP_J,
	'H' = ZP_H,
	'Ks' = ZP_K)

#
# Read in reference files [Planck, GAIA, GAMA IC, GAMA Z, REGROUP, EYEBALL CLASSIFICATIONS]
#
ebv=fread(paste0(RefDir,"ebvtrim.csv"))
gaia=fread(paste0(RefDir,"gaiastarmaskwaves.csv"))

AstroVerification = fread(paste0(RefDir,"gaia_AstroVerification_combined.csv"))

registerDoParallel(cores=CoreNumber)

foreach(j=(1:length(InputTargetCat$ra)))%dopar%{
# for (j in 1:length(InputTargetCat$ra)){
# for (j in c(1)){
	ra=as.character((format(round(InputTargetCat$ra[j], 1), nsmall = 1)))
	dec=as.character((format(round(InputTargetCat$dec[j], 1), nsmall = 1)))
	
	if((as.double(ra) < 12) & (as.double(ra) > 6)){
		region = 'WD01'
	}else if((as.double(ra) < 38) & (as.double(ra) > 33)){
		region = 'WD02'
	}else if((as.double(ra) < 56) & (as.double(ra) > 51)){
		region = 'WD03'
	}else if((as.double(ra) < 152) & (as.double(ra) > 148)){
		region = 'WD10'
	}
	
	
	message(paste0('coordinates: ', ra, ' ', dec))
	
	
	#
	# Read in RDS file to be post-processed
	#
	if(file.exists(paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds')) & !file.exists(paste0(PostprocessDir,"waves_postprocessed_", ra, '_', dec, '.rds'))){
		MeasureFilename=paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds')
		everything=readRDS(MeasureFilename)
		message(paste0('measure file read in: ', ra, ' ', dec))
		#
		# Extract segment info, colour, total, deblend, aperture, and groups measurements
		#
		cat_objects <- as.data.table(cbind(everything$pro_detect$segstats,everything$cat_tot,everything$cat_col)) # including the colour outputs for the sake of photometric redshifts
		cat_groupinfo=cbind(segID=unlist(everything$pro_detect$group$groupsegID$segID), 
			groupID=rep(everything$pro_detect$group$groupsegID$groupID, times=everything$pro_detect$group$groupsegID$Ngroup), 
			Ngroup=rep(everything$pro_detect$group$groupsegID$Ngroup, times=everything$pro_detect$group$groupsegID$Ngroup))
		cat_objects=cbind(cat_objects,cat_groupinfo[match(cat_objects$segID, cat_groupinfo[,"segID"]),2:3])
		cat_groups <- as.data.table(cbind(everything$pro_detect$group$groupsegID$Ngroup,everything$pro_detect$groupstats$groupID,everything$cat_grp))
		names(cat_groups)[1] <- "Ngroup"
		names(cat_groups)[2] <- "groupID"
		group_matches=match(cat_objects$segID,cat_groups$groupID,nomatch=NA)
		datafile0=as.data.table(cbind(cat_objects,cat_groups[group_matches,]))
		message(paste0('coordinates: ', ra, ' ', dec))
		
		message('Add Naming info and estimate frame seeing')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Add Naming info and estimate frame seeing
		#
		framename=paste0("waves_postprocessed_", ra, '_', dec, '.rds')
		# iauname=IAUID(datafile0$RAmax,datafile0$Decmax)
		# frameid=ceiling(min(datafile0$RAmax))*100.0+ceiling(min(datafile0$Decmax)) # ceiling value of ra at min x
		# frameid=ceiling(datafile0$RAmax[which(datafile0$xmax==min(datafile0$xmax))])*100.0+ceiling(min(datafile0$Decmax))
		frameid=rep(ceiling(as.double(ra))*100.0+ceiling(as.double(dec)), length(datafile0$uniqueID))
		frameid=ceiling(as.double(ra))*100.0+ceiling(as.double(dec))
		# uberid = as.integer64(frameid)*1e10 + as.integer64(datafile0$uniqueID)
		uberid=as.character(frameid*1e10+datafile0$uniqueID)
		# print(uberid)
		seeing=log10(median(datafile0[mag_Zt > 16.0 & mag_Zt < 17 & R50 < 2.0,R50]))
		seeing_r=log10(median(datafile0[mag_Zt > 16.0 & mag_Zt < 17 & R50 < 2.0,R50_rt]))
		seeing_i=log10(median(datafile0[mag_Zt > 16.0 & mag_Zt < 17 & R50 < 2.0,R50_it]))
		seeing_Z=log10(median(datafile0[mag_Zt > 16.0 & mag_Zt < 17 & R50 < 2.0,R50_Zt]))
		seeing_Y=log10(median(datafile0[mag_Zt > 16.0 & mag_Zt < 17 & R50 < 2.0,R50_Yt]))
		datafile0$uberID=uberid
		datafile0$FrameName=framename
		datafile0$FrameID=frameid
		datafile0$log10seeing=seeing
		datafile0$log10seeing_r=seeing_r
		datafile0$log10seeing_i=seeing_i
		datafile0$log10seeing_Z=seeing_Z
		datafile0$log10seeing_Y=seeing_Y

		# now correcting the fluxes for the zeropoint offsets
		for(filter in c('u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks')){
			for(photom_type in c('t', 'c')){ 
				zp = as.numeric(ZP_corrections[filter][[1]][region])
				
				band=paste0("mag_",filter,photom_type)
				datafile0[[band]] = datafile0[[band]] + zp
				#
				band=paste0("flux_",filter,photom_type)
				datafile0[[band]] = datafile0[[band]] * (10^(0.4*(-zp)))
				#
				band=paste0("flux_err_",filter,photom_type)
				datafile0[[band]] = datafile0[[band]] * (10^(0.4*(-zp)))
			}
		}

		
		message('Trim Planck and GAIA maps to appropriate region (saves time)')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Trim Planck and GAIA maps to appropriate region (saves time)
		#
		wavesra=c(min(datafile0[,RAmax]),max(datafile0[,RAmax]))
		wavesdec=c(min(datafile0[,Decmax]),max(datafile0[,Decmax]))
		gaia0=gaia[ra > wavesra[1]-1.0 & ra < wavesra[2]+1.0 & dec > wavesdec[1]-1.0 & dec < wavesdec[2]+1.0,]
		AstroVerification0=AstroVerification[ra > wavesra[1]-1.0 & ra < wavesra[2]+1.0 & dec > wavesdec[1]-1.0 & dec < wavesdec[2]+1.0,]
		ebv0=ebv[RA_J2000 > wavesra[1]-1.0 & RA_J2000 < wavesra[2]+1.0 & DEC_J2000 > wavesdec[1]-1.0 & DEC_J2000 < wavesdec[2]+1.0,]
		
		message('Determine extinction (E(B-V) values using Planck E(B-V) map)')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Determine extinction (E(B-V) values using Planck E(B-V) map)
		#
		plank=coordmatch(coordref=datafile0[,list(RAmax,Decmax)],coordcompare=ebv0[,list(RA_J2000,DEC_J2000)],rad=150,inunitref="deg",inunitcompare="deg",radunit="asec")
		datafile0=cbind(datafile0,ebv0[plank$ID[,1],"EBV"])
		
		message('Implement extinction corrections for SBs, mags, fluxes and flux errors')
		message(paste0('coordinates: ', ra, ' ', dec))
		
		# include the non-extinction corrected magnitudes separately in the catalogue, for the sake of target catalogues in 4MOST
		#
		# remove bicubic sky and add local sky
		#  
		# Implement extinction corrections for SBs, mags, fluxes and flux errors
		#
		# Af=8.24152
		# An=8.20733
		Au=4.81139
		Ag=3.66469
		Ar=2.65460
		Ai=2.07472
		AZ=1.55222
		AY=1.21291
		AJ=0.87624
		AH=0.56580
		AK=0.36888
		A1=0.20124
		A2=0.13977
		#
		extinc_coeffs=rbind(Au,Ag,Ar,Ai,AZ,AY,AJ,AH,AK,A1,A2)
		filters=rbind("u","g","r","i","Z","Y","J","H","Ks","W1","W2")
		zeropoints=rbind(23.9,23.9,23.9,23.9,23.9,23.9,23.9,23.9,23.9,23.183,22.819) # updated zeropoint compared to Wide
		#
		bands=as.data.table(cbind(filters,extinc_coeffs,zeropoints))
		names(bands)[1] <- "filters"
		names(bands)[2] <- "extinc_coeffs"
		names(bands)[3] <- "zeropoints"
		#
		for(filter in bands$filters){
			ext=as.numeric(bands[filter==filters,"extinc_coeffs"])
			zp=as.numeric(bands[filter==filters,"zeropoints"])
			#
			sky1=paste0("sky_sum_",filter,"t")
			sky2=paste0("skyseg_mean_",filter,"t")
			Npix=paste0("N100_",filter,"t")
		
			# sky1jy=datafile0[[sky1]]*10^(0.4*(8.9-zp))
			# sky2jy=datafile0[[sky2]]*datafile0[[Npix]]*10^(0.4*(8.9-zp))
			#
			for(photom_type in c('t', 'c')){ # including the colour photometry for photo-z work
				band=paste0("mag_",filter,photom_type)
				datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]-ext*datafile0$EBV)
				#
				band=paste0("mag_app_",filter,photom_type)
				datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]-ext*datafile0$EBV)
				#
				band=paste0("flux_",filter,photom_type)
				# including the uncorrected photometry (added Mar22)
				# needs to be added before corrected fluxes, to avoid overwriting!
				datafile0[[paste0("flux_",filter,photom_type,"_uncorrected")]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]])
				datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))
				#
				# band=paste0("flux_",filter,"l")
				# datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))
				#	
				band=paste0("flux_err_",filter,photom_type)
				# including the uncorrected photometry (added Mar22)
				datafile0[[paste0("flux_err_",filter,photom_type,"_uncorrected")]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]])
				datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))
			}
			
		}
		#
		# Add blank classification columns to datafile
		#
		message('Add blank classification columns to datafile')
		message(paste0('coordinates: ', ra, ' ', dec))

		# detect = Rfits_read(paste0(DetectDir, 'waves_detect_',ra,'_',dec,'.fits'), pointer=FALSE)
		
		
		datafile0 = as.data.table(cbind(datafile0,
			censep=as.numeric(sqrt((12389/2.0-datafile0[,xmax])^2+(12389/2.0-datafile0[,ymax])^2)),
			# RAGAIA=datafile0$RAmax-gaiaraoff/3600.0,
			# DecGAIA=datafile0$Decmax-gaiadecoff/3600.0,
			RAGAIA=datafile0$RAmax,
			DecGAIA=datafile0$Decmax,
			class="notclassified",
			eyeclass="unknown",
			uberclass="unknown",
			# skyRMS_mean_detect = detect$segstats$skyRMS_mean,
			duplicate=1, # not sure if this is the right place to put this - duplicates aren't identified until later. 
			starscol=0,
			starssize=0,
			mask=0,
			noOPT_r=0,
			noOPT_i=0,
			noIR_Z=0,
			noIR_Y=0,
			starmask=0))
		
		message('Assign regions with no optical data')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Assign regions with no optical or IR data
		#
		datafile0[is.na(flux_rt),"noOPT_r"]=1 # need to update what these values actually are for objects with missing data
		datafile0[is.na(flux_it),"noOPT_i"]=1 # need to update what these values actually are for objects with missing data
		datafile0[is.na(flux_Zt),"noIR_Z"]=1 # need to update what these values actually are for objects with missing data
		datafile0[is.na(flux_Yt),"noIR_Y"]=1 # need to update what these values actually are for objects with missing data
		
		message('Assign stellar colour flag ')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Assign stellar colour flag 
		#
		datafile0[mag_rt > 19.5 & (mag_Jt-mag_Kst) < (0.025+0.025*(mag_rt-19.5)),"starscol"]=1 # ambiguous 
		datafile0[mag_rt < 19.5 & ( ((mag_Jt-mag_Kst) < 0.025) | (mag_rt < 12) ),"starscol"]=3 # star - Ah, it looks like this brightness cut is causing some really massive gals to flag as stars!
		datafile0[mag_rt > 19.5 & (mag_Jt-mag_Kst) < (0.025-0.1*(mag_rt-19.5)^2.0),"starscol"]=3 # star
		
		# if there's no r-band data
		datafile0[noOPT_r==1 & (mag_Zt > 18.5) & ((mag_Jt-mag_Kst) < (0.025+0.025*(mag_Zt-18.5))),"starscol"]=1 # ambiguous 
		datafile0[noOPT_r==1 & (mag_Zt < 18.5) & ( ((mag_Jt-mag_Kst) < 0.025) | (mag_Zt < 12) ),"starscol"]=3 # star - Ah, it looks like this brightness cut is causing some really massive gals to flag as stars!
		datafile0[noOPT_r==1 & (mag_Zt > 18.5) & ((mag_Jt-mag_Kst) < (0.025-0.1*(mag_Zt-18.5)^2.0)),"starscol"]=3 # star

		
		
		message('Assign stellar size flag')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		#  Assign stellar size flag
		#
		datafile0[,"starssize"]=0 # anything that stays 0 is in the galaxy section of parameter space
		datafile0[(log10(R50) < seeing+0.05),"starssize"]=1 # ambiguous 		
		datafile0[(log10(R50) < (seeing+0.05-0.075*(mag_Zt-19.5)^1.0)),"starssize"]=3 # star

		# if there's no r-band data
		datafile0[noOPT_r==1 & (log10(R50) < (seeing+0.05-0.075*(mag_Zt-19.5)^1.0)),"starssize"]=3 # star

		
		
		
		message('Assign object classes starting with ambiguous')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Assign object classes starting with ambiguous
		#
		datafile0[,"class"]="ambiguous"
		# now, seting the selection for objects with r, K, and J-band data
		datafile0[(starscol+starssize < 1.5) & !is.na(flux_Jt) & !is.na(flux_Kst),"class"]="galaxy" 
		datafile0[(starscol+starssize > 3.5) & !is.na(flux_Jt) & !is.na(flux_Kst),"class"]="star" 	
		
		# if there's no J or K-band data, then I only want to use the size criterion
		datafile0[(starssize==0) & (is.na(flux_Jt) | is.na(flux_Kst)),"class"]="galaxy" 
		datafile0[(starssize==3) & (is.na(flux_Jt) | is.na(flux_Kst)),"class"]="star" 	
		
		#
		# Identify artefacts (no detection outside r+Z, too small, bright improbable colour, high skyRMS)
		# ah, realistically I probably need to do this in the discosky regions separately...
		#
		# message('Define disco skyRMS region')
		# DiscoSky = autoDisco(skyRMS=everything$pro_detect$skyRMS, blurSigma=200, kernelSmoothing=max(everything$pro_detect$skyRMS)/100)# NEED TO READ IN THIS VALUE HERE
		
		# # now, for each object in the catalogue, identify which sky region its in
		# SkyRegions = DiscoSky[cbind(ceiling(datafile0$xmax), ceiling(datafile0$ymax))]

		# to define the sky regions I actually want to use the original weights maps, because they provide nice sharp edges
		# message('defining sky regions from detection skyRMS (from weights)')
		# DiscoSky = autoDisco(skyRMS=detect$skyRMS$imDat, kernelSmoothing=max(detect$skyRMS$imDat)/100, blur=FALSE)

		# SkyRegions = DiscoSky[cbind(ceiling(datafile0$xmax), ceiling(datafile0$ymax))]

		# # datafile0["skyRMS_mean_detect"]=detect$segstats$skyRMS_mean

		# # saving memory
		# rm(detect)
	
		# PlotFilename=paste0(PlotDir,"postptocess_skyRMSregions_",ra,"_",dec,".png")
		# png(PlotFilename,width=15.0,height=15.0,units="cm",res=120)
		# magimage(DiscoSky)
		# dev.off()
		
		
		# message(paste0('number of skyRMS regions: ', length(unique(SkyRegions))))
		# message(unique(SkyRegions))
		# for(SkyRegion in 1:length(unique(SkyRegions))){
		# 	medianskyRMS=quantile(datafile0[SkyRegions==SkyRegion,skyRMS_mean_detect],0.5)
		# 	lowsky=datafile0[(skyRMS_mean_detect < medianskyRMS) & (SkyRegions==SkyRegion),skyRMS_mean_detect]
		# 	crudcut=medianskyRMS+5*(medianskyRMS-quantile(lowsky,0.33))
			
		# 	message(paste0('	number of artefacts in non-Z data before skyRMS cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR_Z==1])))
		# 	datafile0[(skyRMS_mean_detect > crudcut) & (SkyRegions==SkyRegion),"class"]="artefact" # high sky RMS. Ideally we actually want a measure of the mean skyRMS that just ignores the r band when it's missing...
		# 	# datafile0[skyRMS_mean > crudcut & noOPT_r<1 & (SkyRegions==SkyRegion),"class"]="artefact" # not sure if needed?
		# }

		
		datafile0[is.na(mag_gt)+is.na(mag_rt)+is.na(mag_it) > 2 & noOPT_r < 1,"class"]="artefact" # bright/missing optical, where r-band exists
		# datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kst) > 2 & noOPT_r < 1,"class"]="artefact" # bright/missing NIR, where r-band exists. This finds ghosted regions pretty well. 
		datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kst) > is.na(flux_Zt)+is.na(flux_Yt)+is.na(flux_Jt)+is.na(flux_Ht)+is.na(flux_Kst)+2 & 
							noOPT_r < 1 & noIR_Z<1,"class"]="artefact" # bright/missing NIR, where r-band exists. This finds ghosted regions pretty well. 
		datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kst) > is.na(flux_Zt)+is.na(flux_Yt)+is.na(flux_Jt)+is.na(flux_Ht)+is.na(flux_Kst)+2 & # flux is only na if the data are missing. 
							noOPT_r < 1 & noIR_Z==1,"class"]="artefact" # bright/missing NIR, where r-band exists. This finds ghosted regions pretty well. 
		datafile0[log10(R50) < -0.4,"class"]="artefact" # too small
		datafile0[(mag_rt-mag_Zt) < -0.75 & noIR_Z < 1 & noOPT_r < 1,"class"]="artefact" # improbable r-Z colour
		datafile0[(mag_rt-mag_Zt) > 3 & noIR_Z < 1 & noOPT_r < 1,"class"]="artefact" # improbable r-Z colour (too high)
		datafile0[(mag_it-mag_Zt) < -1.5 & noIR_Z < 1 & noOPT_i < 1,"class"]="artefact" # improbable i-Z colour
		datafile0[(mag_it-mag_Zt) > 3 & noIR_Z < 1 & noOPT_i < 1,"class"]="artefact" # improbable i-Z colour (too high)


		datafile0[(skyRMS_mean) > 0.06,"class"]="artefact" # high skyRMS
		
		# for objects with no r-band data, print out their resulting class:
		message('classes for objects with no r-band data:')
		message(paste0('coordinates: ', ra, ' ', dec))
		message(paste0('number of stars: ', length(datafile0$class[datafile0$class=='star' & datafile0$noOPT_r==1])))
		message(paste0('number of galaxies: ', length(datafile0$class[datafile0$class=='galaxy' & datafile0$noOPT_r==1])))
		message(paste0('number of ambiguous: ', length(datafile0$class[datafile0$class=='ambiguous' & datafile0$noOPT_r==1])))
		message(paste0('number of artefacts: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noOPT_r==1])))

		
		# for objects with no Z-band data, print out their resulting class:
		message('classes for objects with no Z-band data:')
		message(paste0('coordinates: ', ra, ' ', dec))
		message(paste0('number of stars: ', length(datafile0$class[datafile0$class=='star' & datafile0$noIR_Z==1])))
		message(paste0('number of galaxies: ', length(datafile0$class[datafile0$class=='galaxy' & datafile0$noIR_Z==1])))
		message(paste0('number of ambiguous: ', length(datafile0$class[datafile0$class=='ambiguous' & datafile0$noIR_Z==1])))
		message(paste0('number of artefacts: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR_Z==1])))
		
		message('Identify star mask and set starmask flag and GAIA match to star')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Identify star mask and set starmask flag and GAIA match to star
		#
		gaia0=subset(gaia0,gaia0$phot_g_mean_mag<16)
		radius=10^(1.6-0.15*gaia0$phot_g_mean_mag)
		radius=ifelse(gaia0$phot_g_mean_mag<6.0,5.0119,radius)
		#
		maskedregion=coordmatch(coordcompare=datafile0[,list(RAmax,Decmax)],coordref=gaia0[,list(ra,dec)],rad=radius,inunitref="deg",inunitcompare="deg",radunit="amin")
		# maskedregion$bestmatch$refID gives the GAIA row positions
		# maskedregion$bestmatch$compareID gives the WAVES row positions
		# maskedregion$bestmatch$sep gives the separation match in amins. 
		length(maskedregion$bestmatch$refID)
		maskedobjs=unique(maskedregion$ID[maskedregion$ID>0])
		length(maskedobjs)
		datafile0[maskedobjs,"starmask"]=1
		datafile0[maskedregion$bestmatch$compareID,"class"]="star"
		#
		
		# Identify GAMA region for GAIA offsets
		#
		# characterise the RA offsets (added Mar22)
		message('Determine the RA and Dec offsets to GAIA')
		GaiaMatch=coordmatch(coordcompare=datafile0[,list(RAmax,Decmax)],coordref=AstroVerification0[,list(ra,dec)],rad=1.0,inunitref="deg",inunitcompare="deg",radunit="asec")
		message(paste0('number of GAIA matches: ', length(GaiaMatch$bestmatch$compareID)))
		
		gaia_RAoffsets = datafile0$RAcen_rt[GaiaMatch$bestmatch$compareID] - AstroVerification0$ra[GaiaMatch$bestmatch$refID]
		gaia_Decoffsets = datafile0$Deccen_rt[GaiaMatch$bestmatch$compareID] - AstroVerification0$dec[GaiaMatch$bestmatch$refID]
		gaiaraoff=median(gaia_RAoffsets)
		gaiadecoff=median(gaia_Decoffsets)
		message(paste0('median RA offset [deg]: ', gaiaraoff))
		message(paste0('median Dec offset [deg]: ', gaiadecoff))
		
		# now correcting the GAIA coordinates in the catalogue
		datafile0[,"RAGAIA"]=datafile0$RAmax-gaiaraoff
		datafile0[,"DecGAIA"]=datafile0$Decmax-gaiadecoff

		datafile0[,"RAGAIA_r"]=datafile0$RAmax_rt-gaiaraoff
		datafile0[,"DecGAIA_r"]=datafile0$Decmax_rt-gaiadecoff

		datafile0[,"RAGAIA_r_cen"]=datafile0$RAcen_rt-gaiaraoff
		datafile0[,"DecGAIA_r_cen"]=datafile0$Deccen_rt-gaiadecoff
		
		
		# now making a plot that demonstrates the offset to GAIA
		PlotFilename=paste0(PlotDir,"waves_GAIAoffset_",ra,"_",dec,"_DR3.png")
		print(PlotFilename)
		png(PlotFilename,width=15.0,height=15.0,units="cm",res=120)
		magplot(gaia_RAoffsets*3600, gaia_Decoffsets*3600, pch='.', xlab='RA offset [asec]', ylab='Dec offset [asec]', xlim=c(-2, 2), ylim=c(-2, 2))
		abline(v = gaiaraoff*3600, col = "blue", lwd = 1, lty = 4) 
		abline(h = gaiadecoff*3600, col = "blue", lwd = 1, lty = 4) 
		title(main=paste0(region," ",ra," ",dec),
		      sub=paste0('RA offset: ', signif(gaiaraoff*3600, digits = 4), ' +/- ', signif(sd(gaia_RAoffsets*3600), digits = 2), 
		                 '  Dec offset: ', signif(gaiadecoff*3600, digits = 4), ' +/- ', signif(sd(gaia_Decoffsets*3600), digits = 2), 
		                 '  matches: ', length(gaia_RAoffsets)))
		dev.off()

		png(paste0(PlotDir,"waves_GAIAoffset_",ra,"_",dec,"_DR3_multipanel.png"),width=45.0,height=30.0,units="cm",res=120)
		par(mfrow=c(2, 3),mar=c(4,4,4,4))
		gaia_RAoffsets = datafile0$RAmax_gt[GaiaMatch$bestmatch$compareID] - AstroVerification0$ra[GaiaMatch$bestmatch$refID]
		gaia_Decoffsets = datafile0$Decmax_gt[GaiaMatch$bestmatch$compareID] - AstroVerification0$dec[GaiaMatch$bestmatch$refID]
		gaiaraoff=median(gaia_RAoffsets)
		gaiadecoff=median(gaia_Decoffsets)
		magplot(gaia_RAoffsets*3600, gaia_Decoffsets*3600, pch='.', xlab='RA offset [asec]', ylab='Dec offset [asec]', xlim=c(-2, 2), ylim=c(-2, 2))
		abline(v = gaiaraoff*3600, col = "blue", lwd = 1, lty = 4) 
		abline(h = gaiadecoff*3600, col = "blue", lwd = 1, lty = 4) 
		title(main=paste0('max (g) RA offset: ', signif(gaiaraoff*3600, digits = 4), ' +/- ', signif(sd(gaia_RAoffsets*3600), digits = 2), 
		                 '  Dec offset: ', signif(gaiadecoff*3600, digits = 4), ' +/- ', signif(sd(gaia_Decoffsets*3600), digits = 2), 
		                 '  matches: ', length(gaia_RAoffsets)))

		gaia_RAoffsets = datafile0$RAmax_rt[GaiaMatch$bestmatch$compareID] - AstroVerification0$ra[GaiaMatch$bestmatch$refID]
		gaia_Decoffsets = datafile0$Decmax_rt[GaiaMatch$bestmatch$compareID] - AstroVerification0$dec[GaiaMatch$bestmatch$refID]
		gaiaraoff=median(gaia_RAoffsets)
		gaiadecoff=median(gaia_Decoffsets)
		magplot(gaia_RAoffsets*3600, gaia_Decoffsets*3600, pch='.', xlab='RA offset [asec]', ylab='Dec offset [asec]', xlim=c(-2, 2), ylim=c(-2, 2))
		abline(v = gaiaraoff*3600, col = "blue", lwd = 1, lty = 4) 
		abline(h = gaiadecoff*3600, col = "blue", lwd = 1, lty = 4) 
		title(main=paste0('max (r) RA offset: ', signif(gaiaraoff*3600, digits = 4), ' +/- ', signif(sd(gaia_RAoffsets*3600), digits = 2), 
		                 '  Dec offset: ', signif(gaiadecoff*3600, digits = 4), ' +/- ', signif(sd(gaia_Decoffsets*3600), digits = 2), 
		                 '  matches: ', length(gaia_RAoffsets)))

		gaia_RAoffsets = datafile0$RAmax_Zt[GaiaMatch$bestmatch$compareID] - AstroVerification0$ra[GaiaMatch$bestmatch$refID]
		gaia_Decoffsets = datafile0$Decmax_Zt[GaiaMatch$bestmatch$compareID] - AstroVerification0$dec[GaiaMatch$bestmatch$refID]
		gaiaraoff=median(gaia_RAoffsets)
		gaiadecoff=median(gaia_Decoffsets)
		magplot(gaia_RAoffsets*3600, gaia_Decoffsets*3600, pch='.', xlab='RA offset [asec]', ylab='Dec offset [asec]', xlim=c(-2, 2), ylim=c(-2, 2))
		abline(v = gaiaraoff*3600, col = "blue", lwd = 1, lty = 4) 
		abline(h = gaiadecoff*3600, col = "blue", lwd = 1, lty = 4) 
		title(main=paste0('max (Z) RA offset: ', signif(gaiaraoff*3600, digits = 4), ' +/- ', signif(sd(gaia_RAoffsets*3600), digits = 2), 
		                 '  Dec offset: ', signif(gaiadecoff*3600, digits = 4), ' +/- ', signif(sd(gaia_Decoffsets*3600), digits = 2), 
		                 '  matches: ', length(gaia_RAoffsets)))

		gaia_RAoffsets = datafile0$RAcen_gt[GaiaMatch$bestmatch$compareID] - AstroVerification0$ra[GaiaMatch$bestmatch$refID]
		gaia_Decoffsets = datafile0$Deccen_gt[GaiaMatch$bestmatch$compareID] - AstroVerification0$dec[GaiaMatch$bestmatch$refID]
		gaiaraoff=median(gaia_RAoffsets)
		gaiadecoff=median(gaia_Decoffsets)
		magplot(gaia_RAoffsets*3600, gaia_Decoffsets*3600, pch='.', xlab='RA offset [asec]', ylab='Dec offset [asec]', xlim=c(-2, 2), ylim=c(-2, 2))
		abline(v = gaiaraoff*3600, col = "blue", lwd = 1, lty = 4) 
		abline(h = gaiadecoff*3600, col = "blue", lwd = 1, lty = 4) 
		title(main=paste0('cen (g) RA offset: ', signif(gaiaraoff*3600, digits = 4), ' +/- ', signif(sd(gaia_RAoffsets*3600), digits = 2), 
		                 '  Dec offset: ', signif(gaiadecoff*3600, digits = 4), ' +/- ', signif(sd(gaia_Decoffsets*3600), digits = 2), 
		                 '  matches: ', length(gaia_RAoffsets)))

		gaia_RAoffsets = datafile0$RAcen_rt[GaiaMatch$bestmatch$compareID] - AstroVerification0$ra[GaiaMatch$bestmatch$refID]
		gaia_Decoffsets = datafile0$Deccen_rt[GaiaMatch$bestmatch$compareID] - AstroVerification0$dec[GaiaMatch$bestmatch$refID]
		gaiaraoff=median(gaia_RAoffsets)
		gaiadecoff=median(gaia_Decoffsets)
		magplot(gaia_RAoffsets*3600, gaia_Decoffsets*3600, pch='.', xlab='RA offset [asec]', ylab='Dec offset [asec]', xlim=c(-2, 2), ylim=c(-2, 2))
		abline(v = gaiaraoff*3600, col = "blue", lwd = 1, lty = 4) 
		abline(h = gaiadecoff*3600, col = "blue", lwd = 1, lty = 4) 
		title(main=paste0('cen (r) RA offset: ', signif(gaiaraoff*3600, digits = 4), ' +/- ', signif(sd(gaia_RAoffsets*3600), digits = 2), 
		                 '  Dec offset: ', signif(gaiadecoff*3600, digits = 4), ' +/- ', signif(sd(gaia_Decoffsets*3600), digits = 2), 
		                 '  matches: ', length(gaia_RAoffsets)))

		gaia_RAoffsets = datafile0$RAcen_Zt[GaiaMatch$bestmatch$compareID] - AstroVerification0$ra[GaiaMatch$bestmatch$refID]
		gaia_Decoffsets = datafile0$Deccen_Zt[GaiaMatch$bestmatch$compareID] - AstroVerification0$dec[GaiaMatch$bestmatch$refID]
		gaiaraoff=median(gaia_RAoffsets)
		gaiadecoff=median(gaia_Decoffsets)
		magplot(gaia_RAoffsets*3600, gaia_Decoffsets*3600, pch='.', xlab='RA offset [asec]', ylab='Dec offset [asec]', xlim=c(-2, 2), ylim=c(-2, 2))
		abline(v = gaiaraoff*3600, col = "blue", lwd = 1, lty = 4) 
		abline(h = gaiadecoff*3600, col = "blue", lwd = 1, lty = 4) 
		title(main=paste0('cen (Z) RA offset: ', signif(gaiaraoff*3600, digits = 4), ' +/- ', signif(sd(gaia_RAoffsets*3600), digits = 2), 
		                 '  Dec offset: ', signif(gaiadecoff*3600, digits = 4), ' +/- ', signif(sd(gaia_Decoffsets*3600), digits = 2), 
		                 '  matches: ', length(gaia_RAoffsets)))
		dev.off()
		
		
		
		# set all objects within 50 pixels of frame edge as artefact
		#
		#datafile0[xmax < 850 | xmax > 11488 | ymax < 850 | ymax > 11488,"class"]="atedge"
		#
		# renaming the bands with Ks to match the wide format
		colnames(datafile0)[colnames(datafile0) == "flux_Kst"] <- "flux_Kt"
		colnames(datafile0)[colnames(datafile0) == "flux_err_Kst"] <- "flux_err_Kt"
		colnames(datafile0)[colnames(datafile0) == "flux_Kst_uncorrected"] <- "flux_Kt_uncorrected"
		colnames(datafile0)[colnames(datafile0) == "flux_err_Kst_uncorrected"] <- "flux_err_Kt_uncorrected"
		colnames(datafile0)[colnames(datafile0) == "flux_Ksc"] <- "flux_Kc"
		colnames(datafile0)[colnames(datafile0) == "flux_err_Ksc"] <- "flux_err_Kc"
		
		#
		# Only output final columns to save space
		#
		fdatafile0=datafile0[,
		  c("FrameName","FrameID","uberID","segID","xmax","ymax","censep",
		  	"RAcen","Deccen","RAmax","Decmax","RAGAIA_r","RAGAIA","DecGAIA","DecGAIA_r","RAGAIA_r_cen","DecGAIA_r_cen",
		  	"RAcen_gt","Deccen_gt","RAmax_gt","Decmax_gt",
		  	"RAcen_rt","Deccen_rt","RAmax_rt","Decmax_rt",
		  	"RAcen_Zt","Deccen_Zt","RAmax_Zt","Decmax_Zt",
		    "sky_mean","skyRMS_mean",
		    "log10seeing","log10seeing_r", "log10seeing_i", "log10seeing_Z", "log10seeing_Y",
		    "mag","EBV",
		    "R50","R50_gt","R50_rt","R50_Yt","R50_Jt","R50_Ht",
		    "R90","R100","N100","axrat","ang",
		    "groupID","Ngroup","mag_app_Zt","mag_Zt",
		    "mag_app_gt","mag_app_rt","mag_app_it","mag_app_Yt",
		    # "flux_FUVt","flux_FUVl","flux_err_FUVt",
		    # "flux_NUVt","flux_NUVl","flux_err_NUVt",
		    "flux_ut","flux_err_ut","flux_ut_uncorrected","flux_err_ut_uncorrected",#"flux_ul",
		    "flux_gt","flux_err_gt","flux_gt_uncorrected","flux_err_gt_uncorrected",#"flux_gl",
		    "flux_rt","flux_err_rt","flux_rt_uncorrected","flux_err_rt_uncorrected",#"flux_rl",
		    "flux_it","flux_err_it","flux_it_uncorrected","flux_err_it_uncorrected",#"flux_il",
		    "flux_Zt","flux_err_Zt","flux_Zt_uncorrected","flux_err_Zt_uncorrected",#"flux_Zl",
		    "flux_Yt","flux_err_Yt","flux_Yt_uncorrected","flux_err_Yt_uncorrected",#"flux_Yl",
		    "flux_Jt","flux_err_Jt","flux_Jt_uncorrected","flux_err_Jt_uncorrected",#"flux_Jl",
		    "flux_Ht","flux_err_Ht","flux_Ht_uncorrected","flux_err_Ht_uncorrected",#"flux_Hl",
		    "flux_Kt","flux_err_Kt","flux_Kt_uncorrected","flux_err_Kt_uncorrected",#"flux_Kl",
		    "flux_W1t","flux_err_W1t","flux_W1t_uncorrected","flux_err_W1t_uncorrected",#"flux_W1l",
		    "flux_W2t","flux_err_W2t","flux_W2t_uncorrected","flux_err_W2t_uncorrected",#"flux_W2l",
		    "flux_uc","flux_err_uc",
		    "flux_gc","flux_err_gc",
		    "flux_rc","flux_err_rc",
		    "flux_ic","flux_err_ic",
		    "flux_Zc","flux_err_Zc",
		    "flux_Yc","flux_err_Yc",
		    "flux_Jc","flux_err_Jc",
		    "flux_Hc","flux_err_Hc",
		    "flux_Kc","flux_err_Kc",
		    "flux_W1c","flux_err_W1c",
		    "flux_W2c","flux_err_W2c",
		    # "flux_W1t","flux_W1l","flux_err_W1t",
		    # "flux_W2t","flux_W2l","flux_err_W2t",
		    "mask","starmask","starscol","starssize","class",
		    "noOPT_r", "noOPT_i","noIR_Z", "noIR_Y")]
		
		message('Construct final data structure')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Construct final data sructure
		#
		trim=list()
		trim$cat <- fdatafile0
		trim$detect_segim=everything$pro_detect$detect_segim
		trim$dilated_segim=everything$pro_detect$dilated_segim
		trim$fixed_segim=everything$pro_detect$segim_orig
		trim$fixed_dilated_segim=everything$pro_detect$segim
		trim$group_segim=everything$pro_detect$group$groupim
		trim$final_segim=everything$segimlist
		trim$header=everything$pro_detect$header
		trim$segID_merge=everything$pro_detect$segID_merge
		
		message('save trimmed RDS file')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# save trimmed RDS file
		#
		saveRDS(trim,file=paste0(PostprocessDir,framename))
		#
		# PlotFilename=paste0(PlotDir,"waves_objects_",ra,"_",dec,"_WISE.png")
		# print(PlotFilename)
		# png(PlotFilename,width=25.0,height=25.0,units="cm",res=240)
		# #
		# Rwcs_image(everything$pro_detect$image, header=everything$pro_detect$header)
		# points(fdatafile0[starmask > 0,list(xmax,ymax)],pch=16,col="red",cex=0.25)
		# # points(fdatafile0[CATAID > 0 & Z > 0.002 & NQ > 2 & starmask < 1,list(xmax,ymax)],pch=".",col="green",cex=0.25)
		# points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & class=="artefact",list(xmax,ymax)],pch=4,col="cyan",cex=0.25)
		# # points(fdatafile0[noIR_Z>0,list(xmax,ymax)],pch=18,col="brown4",cex=0.5)
		# # points(fdatafile0[noOPT_r>0,list(xmax,ymax)],pch=18,col="blue4",cex=0.5)
		# # points(fdatafile0[noIR_Y>0,list(xmax,ymax)],pch=18,col="chocolate4",cex=0.5)
		# # points(fdatafile0[noOPT_i>0,list(xmax,ymax)],pch=18,col="aquamarine4",cex=0.5)
		# points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) &  class=="star",list(xmax,ymax)],pch=18,col="yellow",cex=0.5)
		# # points(fdatafile0[flux_Zt < 10^(0.4*(8.9-21.1)) & starmask < 1 & class!="star",list(xmax,ymax)],pch=20,col="gray",cex=0.25)
		# points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="ambiguous",list(xmax,ymax)],pch=1,col="gray",cex=0.75)
		# points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="galaxy",list(xmax,ymax)],pch=5,col="purple",cex=0.75)
		# points(fdatafile0[flux_it > 10^(0.4*(8.9-22)) & starmask < 1 & class=="galaxy" & noIR_Z==1,list(xmax,ymax)],pch=5,col="blue",cex=0.75) # if Z-band data are missing, plot the i<22 gals
		# #
		# gaia0=subset(gaia0,gaia0$phot_g_mean_mag<16)
		# radius=10^(1.6-0.15*gaia0$phot_g_mean_mag)
		# radius=ifelse(gaia0$phot_g_mean_mag<6.0,5.0119,radius)
		# radius=radius*60.0/0.339
		# #
		# xy=Rwcs_s2p(RA=gaia0$ra,Dec=gaia0$dec,header=everything$pro_detect$header)
		# #
		# for (i in 1:length(radius)){
		#   x=xy[i,1]
		#   y=xy[i,2]
		#   r=radius[i]
		#   draw.circle(x,y,r,border="yellow",lty=1,lwd=1)
		# #  cat(x,y,r,"\n")
		# }
		# message("Z<21.1 gals: ",length(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="galaxy",segID]))
		# message("Artefacts: ",length(fdatafile0[class=="artefact",segID]))
		# message("i<22 gals (no Z): ",length(fdatafile0[flux_it > 10^(0.4*(8.9-22)) & starmask < 1 & class=="galaxy" & noIR_Z==1,segID]))
		# message("Masked: ",length(fdatafile0[starmask > 0,segID]))
		# message("noOPT_r: ",length(fdatafile0[noOPT_r > 0,segID]))
		# message("noIR_Z: ",length(fdatafile0[noIR_Z > 0,segID]))
		# title(main=paste0(region," ",ra," ",dec),
		#       sub=paste0("Z<21.1 gals:",length(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="galaxy",segID]),
		#       					 # "     i<22 gals (no Z):",length(fdatafile0[flux_it > 10^(0.4*(8.9-22)) & starmask < 1 & class=="galaxy" & noIR_Z==1,segID]),
		#                  "     Artefacts: ",length(fdatafile0[class=="artefact",segID]),
		#                  "     Masked: ",length(fdatafile0[starmask > 0,segID]),
		#                  "     noOPT_r: ",length(fdatafile0[noOPT_r > 0,segID]),
		#                  "     noOPT_i: ",length(fdatafile0[noOPT_i > 0,segID]),
		#                  "     noIR_Z: ",length(fdatafile0[noIR_Z > 0,segID]),
		#                  "     noIR_Y: ",length(fdatafile0[noIR_Y > 0,segID])))
		# #
		# dev.off()




		#
		
		###
	}else{
		message('missing measure output')
	}
}	