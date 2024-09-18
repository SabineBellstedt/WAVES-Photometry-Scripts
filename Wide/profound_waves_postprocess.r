# .libPaths(c('/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library/',.libPaths()))
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
CoreNumber = 4


#
#  DEFINE FUNCTIONs
#
magAB2Jansky=function(x){10^(-0.4*(x-8.9))} # taken from ProSpect
#


################
#
# MAIN CODE
#
################

region = 'WAVES-S'
# region = 'G09'
# region = 'WAVES-N'
# region = 'test'


WAVESwideDir = '/Volumes/ThunderBay/WAVES/'

MeasureDir = paste0(WAVESwideDir, "profound/measure_0.3Res_WISE/")
PostprocessDir = paste0(WAVESwideDir, "profound/postprocess_0.3Res_WISE/")
PlotDir = paste0(WAVESwideDir, "profound/plots/")
# FixesDir = paste0(WAVESwideDir, "profound/segmentFixes/")

RefDir=paste0(WAVESwideDir, "ref_Sabine/")

InputTargetCat=fread(paste0(RefDir, 'Tiles_', region, '.txt'))

#
# Read in reference files [Planck, GAIA, GAMA IC, GAMA Z, REGROUP, EYEBALL CLASSIFICATIONS]
#

ebv=fread(paste0(RefDir,"ebvtrim.csv"))
gaia=fread(paste0(RefDir,"gaiastarmaskwaves.csv"))

AstroVerification = fread(paste0(RefDir,"gaia_AstroVerification_combined.csv"))

registerDoParallel(cores=CoreNumber)

foreach(j=(1:length(InputTargetCat$RA)))%dopar%{
# for(j in c(1)){

	ra=as.character((format(round(InputTargetCat$RA[j], 1), nsmall = 1)))
	dec=as.character((format(round(InputTargetCat$Dec[j], 1), nsmall = 1)))

	message(paste0('coordinates: ', ra, ' ', dec))
	#
	# Read in RDS file to be post-processed
	#
	if(file.exists(paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds'))  & (!file.info(paste0(PostprocessDir,"waves_postprocessed_", ra, '_', dec, '.rds'))$mtime > "2024-04-05 00:00:00 AWST")){
		MeasureFilename=paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds')
		everything=readRDS(MeasureFilename)
		#
		# Extract segment info, colour, total, deblend, aperture, and groups measurements
		#
		cat_objects <- as.data.table(cbind(everything$pro_detect$segstats,everything$cat_tot,everything$cat_col))
		cat_groupinfo=cbind(segID=unlist(everything$pro_detect$group$groupsegID$segID), 
			groupID=rep(everything$pro_detect$group$groupsegID$groupID, times=everything$pro_detect$group$groupsegID$Ngroup), 
			Ngroup=rep(everything$pro_detect$group$groupsegID$Ngroup, times=everything$pro_detect$group$groupsegID$Ngroup))
		cat_objects=cbind(cat_objects,cat_groupinfo[match(cat_objects$segID, cat_groupinfo[,"segID"]),2:3])
		cat_groups <- as.data.table(cbind(everything$pro_detect$group$groupsegID$Ngroup,everything$pro_detect$groupstats$groupID,everything$cat_grp))
		names(cat_groups)[1] <- "Ngroup"
		names(cat_groups)[2] <- "groupID"
		group_matches=match(cat_objects$segID,cat_groups$groupID,nomatch=NA)
		datafile0=as.data.table(cbind(cat_objects,cat_groups[group_matches,]))
		
		message('Add Naming info and estimate frame seeing')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Add Naming info and estimate frame seeing
		#
		framename=paste0("waves_postprocessed_", ra, '_', dec, '.rds')
		# iauname=IAUID(datafile0$RAmax,datafile0$Decmax)
		# frameid=ceiling(min(datafile0$RAmax))*100.0+ceiling(min(datafile0$Decmax)) # ceiling value of ra at min x
		frameid=ceiling(datafile0$RAmax[which(datafile0$xmax==min(datafile0$xmax))][1])*100.0+ceiling(min(datafile0$Decmax))
		# frameid=rep(ceiling(as.double(ra))*100.0+ceiling(as.double(dec)), length(datafile0$uniqueID))
		# frameid=ceiling(as.double(ra))*100.0+ceiling(as.double(dec)) # this has issues for -ve RA values, whereas RAmax goes to high 300s...
		print(frameid)
		print(ceiling(as.double(ra))*100.0+ceiling(as.double(dec)))
		uberid=as.character(frameid*1e10+datafile0$uniqueID)
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
		filters=rbind("u","g","r","i","Z","Y","J","H","K","W1","W2")
		zeropoints=rbind(0,0,0,0,30,30,30,30,30,23.183,22.819)
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
			band1=paste0("flux_",filter,"t")
			# band2=paste0("flux_",filter,"l")
			sky1=paste0("sky_sum_",filter,"t")
			sky2=paste0("skyseg_mean_",filter,"t")
			Npix=paste0("N100_",filter,"t")
			sky1jy=datafile0[[sky1]]*10^(0.4*(8.9-zp))
			sky2jy=datafile0[[sky2]]*datafile0[[Npix]]*10^(0.4*(8.9-zp))
			# datafile0[[band2]]=ifelse(datafile0[[band1]]==-999,-999,datafile0[[band1]]+sky1jy-sky2jy)
			 # datafile0$flux_Kt + datafile0$sky_sum_Kt - datafile0$skyseg_mean_Kt*datafile0$N100_Kt
			#
			band=paste0("mag_",filter,"t")
			datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]-ext*datafile0$EBV)
			#
			band=paste0("flux_",filter,"t")
			# including the uncorrected photometry (added Mar22)
			# needs to be added before corrected fluxes, to avoid overwriting!
			datafile0[[paste0("flux_",filter,"t_uncorrected")]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]])
			datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))
			#
			# band=paste0("flux_",filter,"l")
			# datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))
			#	
			band=paste0("flux_err_",filter,"t")
			# including the uncorrected photometry (added Mar22)
			datafile0[[paste0("flux_err_",filter,"t_uncorrected")]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]])
			datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))

			# now adding the colour photometry
			band=paste0("mag_",filter,"c")
			datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]-ext*datafile0$EBV)
			#
			band=paste0("flux_",filter,"c")
			datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))
			#	
			band=paste0("flux_err_",filter,"c")
			datafile0[[band]]=ifelse(datafile0[[band]]==-999,-999,datafile0[[band]]*10^(0.4*(ext*datafile0$EBV)))
		}
		#
		# Add blank classification columns to datafile
		#
		message('Add blank classification columns to datafile')
		message(paste0('coordinates: ', ra, ' ', dec))
		
		datafile0 = as.data.table(cbind(datafile0,
			censep=as.numeric(sqrt((12389/2.0-datafile0[,xmax])^2+(12389/2.0-datafile0[,ymax])^2)),
			# RAGAIA=datafile0$RAmax-gaiaraoff/3600.0,
			# DecGAIA=datafile0$Decmax-gaiadecoff/3600.0,
			RAGAIA=datafile0$RAmax,
			DecGAIA=datafile0$Decmax,
			class="notclassified",
			duplicate=1, # not sure if this is the right place to put this - duplicates aren't identified until later. 
			starscol=0,
			starssize=0,
			mask=0,
			noOPT=0,
			noIR=0,
			starmask=0))
		
		message('Assign regions with no optical data')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
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
		datafile0[mag_rt > 19.5 & (mag_Jt-mag_Kt) < (0.025+0.025*(mag_rt-19.5)),"starscol"]=1 # ambiguous 
		datafile0[mag_rt < 19.5 & ( ((mag_Jt-mag_Kt) < 0.025) | (mag_rt < 12) ),"starscol"]=3 # star - Ah, it looks like this brightness cut is causing some really massive gals to flag as stars!
		datafile0[mag_rt > 19.5 & (mag_Jt-mag_Kt) < (0.025-0.1*(mag_rt-19.5)^2.0),"starscol"]=3 # star
		
		# if there's no r-band data
		datafile0[noOPT_r==1 & (mag_Zt > 18.5) & ((mag_Jt-mag_Kt) < (0.025+0.025*(mag_Zt-18.5))),"starscol"]=1 # ambiguous 
		datafile0[noOPT_r==1 & (mag_Zt < 18.5) & ( ((mag_Jt-mag_Kt) < 0.025) | (mag_Zt < 12) ),"starscol"]=3 # star - Ah, it looks like this brightness cut is causing some really massive gals to flag as stars!
		datafile0[noOPT_r==1 & (mag_Zt > 18.5) & ((mag_Jt-mag_Kt) < (0.025-0.1*(mag_Zt-18.5)^2.0)),"starscol"]=3 # star
		
		
		message('Assign stellar size flag')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		#  Assign stellar size flag
		#
		datafile0[,"starssize"]=0 # anything that stays 0 is in the galaxy section of parameter space
		datafile0[(log10(R50) < seeing+0.05),"starssize"]=1 # ambiguous 
		datafile0[(log10(R50) < (seeing+0.05-0.075*(mag_rt-20.5)^1.0)),"starssize"]=3 # star
		
		# if there's no r-band data
		datafile0[noOPT==1 & (log10(R50) < (seeing+0.05-0.075*(mag_Zt-19.5)^1.0)),"starssize"]=3 # star
		
		
		
		message('Assign object classes starting with ambiguous')
		message(paste0('coordinates: ', ra, ' ', dec))
		#
		# Assign object classes starting with ambiguous
		#
		datafile0[,"class"]="ambiguous"
		# now, seting the selection for objects with r, K, and J-band data
		datafile0[(starscol+starssize < 1.5) & !is.na(flux_Jt) & !is.na(flux_Kt),"class"]="galaxy" 
		datafile0[(starscol+starssize > 3.5) & !is.na(flux_Jt) & !is.na(flux_Kt),"class"]="star" 	
		
		# if there's no J or K-band data, then I only want to use the size criterion
		datafile0[(starssize==0) & (is.na(flux_Jt) | is.na(flux_Kt)),"class"]="galaxy" 
		datafile0[(starssize==3) & (is.na(flux_Jt) | is.na(flux_Kt)),"class"]="star" 	
		
		#
		# Identify artefacts (no detection outside r+Z, too small, bright improbable colour, high skyRMS)
		#
		medianskyRMS=quantile(datafile0[,skyRMS_mean],0.5)
		lowsky=datafile0[skyRMS_mean < medianskyRMS,skyRMS_mean]
		crudcut=medianskyRMS+5*(medianskyRMS-quantile(lowsky,0.33))
		
		message(paste0('	number of artefacts in non-r data before skyRMS cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noOPT==1])))
		message(paste0('	number of artefacts in non-Z data before skyRMS cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR==1])))
		# datafile0[skyRMS_mean > crudcut,"class"]="artefact" # high sky RMS
		datafile0[skyRMS_mean > crudcut & noOPT<1,"class"]="artefact" # high sky RMS. Ideally we actually want a measure of the mean skyRMS that just ignores the r band when it's missing...
		message(paste0('	number of artefacts in non-r data before bright/missing optical cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noOPT==1])))
		message(paste0('	number of artefacts in non-Z data before bright/missing optical cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR==1])))
		datafile0[is.na(mag_gt)+is.na(mag_rt)+is.na(mag_it) > 2 & noOPT < 1,"class"]="artefact" # bright/missing optical, where r-band exists
		message(paste0('	number of artefacts in non-r data before bright/missing NIR cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noOPT==1])))
		message(paste0('	number of artefacts in non-Z data before bright/missing NIR cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR==1])))
		# datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kt) > 2 & noOPT < 1,"class"]="artefact" # bright/missing NIR, where r-band exists. This finds ghosted regions pretty well. 
		datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kt) > is.na(flux_Zt)+is.na(flux_Yt)+is.na(flux_Jt)+is.na(flux_Ht)+is.na(flux_Kt)+2 & 
							noOPT < 1 & noIR<1,"class"]="artefact" # bright/missing NIR, where r-band exists. This finds ghosted regions pretty well. 
		datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kt) > is.na(flux_Zt)+is.na(flux_Yt)+is.na(flux_Jt)+is.na(flux_Ht)+is.na(flux_Kt)+2 & # flux is only na if the data are missing. 
							noOPT < 1 & 
							noIR==1,"class"]="artefact" # bright/missing NIR, where r-band exists. This finds ghosted regions pretty well. 
		message(paste0('	number of artefacts in non-r data before size cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noOPT==1])))
		message(paste0('	number of artefacts in non-Z data before size cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR==1])))
		datafile0[log10(R50) < -0.4,"class"]="artefact" # too small
		message(paste0('	number of artefacts in non-r data before colour cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noOPT==1])))
		message(paste0('	number of artefacts in non-Z data before colour cut: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR==1])))
		datafile0[(mag_rt-mag_Zt) < -0.75 & noIR < 1 & noOPT < 1,"class"]="artefact" # improbably r-Z colour
		
		# for objects with no r-band data, print out their resulting class:
		message('classes for objects with no r-band data:')
		message(paste0('coordinates: ', ra, ' ', dec))
		message(paste0('number of stars: ', length(datafile0$class[datafile0$class=='star' & datafile0$noOPT==1])))
		message(paste0('number of galaxies: ', length(datafile0$class[datafile0$class=='galaxy' & datafile0$noOPT==1])))
		message(paste0('number of ambiguous: ', length(datafile0$class[datafile0$class=='ambiguous' & datafile0$noOPT==1])))
		message(paste0('number of artefacts: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noOPT==1])))
		
		# for objects with no Z-band data, print out their resulting class:
		message('classes for objects with no Z-band data:')
		message(paste0('coordinates: ', ra, ' ', dec))
		message(paste0('number of stars: ', length(datafile0$class[datafile0$class=='star' & datafile0$noIR==1])))
		message(paste0('number of galaxies: ', length(datafile0$class[datafile0$class=='galaxy' & datafile0$noIR==1])))
		message(paste0('number of ambiguous: ', length(datafile0$class[datafile0$class=='ambiguous' & datafile0$noIR==1])))
		message(paste0('number of artefacts: ', length(datafile0$class[datafile0$class=='artefact' & datafile0$noIR==1])))
		
		
		# message('Override classification if quality redshift')
		# #
		# # Override classification if quality redshift
		# #
		# datafile0$class=ifelse(datafile0$Z>0.002 & datafile0$NQ>2,"galaxy",datafile0$class)
		# datafile0$class=ifelse(datafile0$Z>-0.002 & datafile0$Z<0.002 & datafile0$NQ>2,"star",datafile0$class)
		
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
		
		# now masking the regions around the two bright globular clusters
		# first masking NGC 1049 at ra=39.925, -34.2833, with a masking radius of 0.5 deg
		# GCcoords= as.data.frame(list(ra=c(39.925, 15.0), dec=c(-34.2833, -33.716)))
		# maskedregion=coordmatch(coordref=GCcoords,coordcompare=datafile0[,list(RAmax,Decmax)],rad=c(0.5, 0.4),inunitref="deg",inunitcompare="deg",radunit="deg")
		# message(paste0('objects masked for NGC 1049 and other GC: ', length(datafile0$RAmax[maskedregion$bestmatch$compareID]) ))
		# print(datafile0$RAmax[maskedregion$bestmatch$compareID])
		# maskedobjs=unique(maskedregion$ID[maskedregion$ID>0])
		# datafile0[maskedobjs,"mask"]=1 # actively masking the region
		#
		# masking the GCs individually
		maskedregion=coordmatchsing(RAref=39.925,Decref=-34.6,coordcompare=datafile0[,list(RAmax,Decmax)],rad=0.7,inunitref="deg",inunitcompare="deg",radunit="deg")
		maskedobjs=unique(maskedregion$ID[maskedregion$ID>0])
		message(paste0('objects masked for NGC 1049: ', length(datafile0$RAmax[maskedobjs]) ))
		datafile0[maskedobjs,"mask"]=1 # actively masking the region
		
		
		# masking the GCs individually
		maskedregion=coordmatchsing(RAref=15.0,Decref=-33.716,coordcompare=datafile0[,list(RAmax,Decmax)],rad=0.35,inunitref="deg",inunitcompare="deg",radunit="deg")
		maskedobjs=unique(maskedregion$ID[maskedregion$ID>0])
		message(paste0('objects masked for other GC: ', length(datafile0$RAmax[maskedobjs]) ))
		datafile0[maskedobjs,"mask"]=1 # actively masking the region
		
		message('Only output final columns to save space')
		message(paste0('coordinates: ', ra, ' ', dec))

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
		# Construct final data structure
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
		PlotFilename=paste0(PlotDir,"waves_objects_",ra,"_",dec,"_WISE.png")
		print(PlotFilename)
		png(PlotFilename,width=25.0,height=25.0,units="cm",res=240)
		#
		Rwcs_image(everything$pro_detect$image, header=everything$pro_detect$header)
		points(fdatafile0[starmask > 0,list(xmax,ymax)],pch=16,col="red",cex=0.25)
		# points(fdatafile0[CATAID > 0 & Z > 0.002 & NQ > 2 & starmask < 1,list(xmax,ymax)],pch=".",col="green",cex=0.25)
		points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & class=="artefact",list(xmax,ymax)],pch=4,col="cyan",cex=0.25)
		# points(fdatafile0[noIR_Z>0,list(xmax,ymax)],pch=18,col="brown4",cex=0.5)
		# points(fdatafile0[noOPT_r>0,list(xmax,ymax)],pch=18,col="blue4",cex=0.5)
		# points(fdatafile0[noIR_Y>0,list(xmax,ymax)],pch=18,col="chocolate4",cex=0.5)
		# points(fdatafile0[noOPT_i>0,list(xmax,ymax)],pch=18,col="aquamarine4",cex=0.5)
		points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) &  class=="star",list(xmax,ymax)],pch=18,col="yellow",cex=0.5)
		# points(fdatafile0[flux_Zt < 10^(0.4*(8.9-21.1)) & starmask < 1 & class!="star",list(xmax,ymax)],pch=20,col="gray",cex=0.25)
		points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="ambiguous",list(xmax,ymax)],pch=1,col="gray",cex=0.75)
		points(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="galaxy",list(xmax,ymax)],pch=5,col="purple",cex=0.75)
		points(fdatafile0[flux_it > 10^(0.4*(8.9-22)) & starmask < 1 & class=="galaxy" & noIR_Z==1,list(xmax,ymax)],pch=5,col="blue",cex=0.75) # if Z-band data are missing, plot the i<22 gals
		#
		gaia0=subset(gaia0,gaia0$phot_g_mean_mag<16)
		radius=10^(1.6-0.15*gaia0$phot_g_mean_mag)
		radius=ifelse(gaia0$phot_g_mean_mag<6.0,5.0119,radius)
		radius=radius*60.0/0.339
		#
		xy=Rwcs_s2p(RA=gaia0$ra,Dec=gaia0$dec,header=everything$pro_detect$header)
		#
		for (i in 1:length(radius)){
		  x=xy[i,1]
		  y=xy[i,2]
		  r=radius[i]
		  draw.circle(x,y,r,border="yellow",lty=1,lwd=1)
		#  cat(x,y,r,"\n")
		}
		message("Z<21.1 gals: ",length(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="galaxy",segID]))
		message("Artefacts: ",length(fdatafile0[class=="artefact",segID]))
		message("i<22 gals (no Z): ",length(fdatafile0[flux_it > 10^(0.4*(8.9-22)) & starmask < 1 & class=="galaxy" & noIR_Z==1,segID]))
		message("Masked: ",length(fdatafile0[starmask > 0,segID]))
		message("noOPT_r: ",length(fdatafile0[noOPT_r > 0,segID]))
		message("noIR_Z: ",length(fdatafile0[noIR_Z > 0,segID]))
		title(main=paste0(region," ",ra," ",dec),
		      sub=paste0("Z<21.1 gals:",length(fdatafile0[flux_Zt > 10^(0.4*(8.9-21.1)) & starmask < 1 & class=="galaxy",segID]),
		      					 # "     i<22 gals (no Z):",length(fdatafile0[flux_it > 10^(0.4*(8.9-22)) & starmask < 1 & class=="galaxy" & noIR_Z==1,segID]),
		                 "     Artefacts: ",length(fdatafile0[class=="artefact",segID]),
		                 "     Masked: ",length(fdatafile0[starmask > 0,segID]),
		                 "     noOPT_r: ",length(fdatafile0[noOPT_r > 0,segID]),
		                 "     noOPT_i: ",length(fdatafile0[noOPT_i > 0,segID]),
		                 "     noIR_Z: ",length(fdatafile0[noIR_Z > 0,segID]),
		                 "     noIR_Y: ",length(fdatafile0[noIR_Y > 0,segID])))
		#
		dev.off()
		#

		rm(everything)
		rm(datafile0)
		rm(fdatafile0)
	}else{
		message(paste0(PostprocessDir,"waves_postprocessed_", ra, '_', dec, '.rds: output exists'))
	}
}
###
