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
library(arrow)
library(dplyr)
library(Rfits)
library(Rwcs)
library(foreach)
library(celestial)
library(ProFound)
#
# Define functions
#
AngleConversion=function(value){
  if (value<0){
    return (360+value)
  }else if(value > 360){
    return (value-360)
  }else{
    return (value)
  }
}
TileNumberCheck=function(RA_coordinates, Dec_coordinates, Frames, FieldRAmin, FieldRAmax, FieldDecmin, FieldDecmax, 
                         ImageDir="/group/pawsey0160/waves/profound/tiles_0.3Res/", start_time){
  # reading all of the relevant headers once, so that I don't need to do it many times
  ImageHeader <- list()
  ImageDimension <- list()
  for(ii in 1:length(Frames)){
    temp=Rfits_read_header(paste(ImageDir, 'r_', Frames[ii], '.fits', sep=''))
    ImageHeader[[ii]]=temp$hdr
    ImageDimension[[ii]]=temp$keyvalues$NAXIS1
  }
  
  TileNumber=rep(0, length(RA_coordinates))
  for(jj in 1:length(RA_coordinates)){
    if(jj%%100000==0){
      print(paste(100*jj/length(RA_coordinates), '%, time elapsed: ', (difftime(Sys.time(), start_time, units = "mins")), ' minutes', sep=''))
    }
    Coordinates=c(RA_coordinates[jj], Dec_coordinates[jj])
    
    # first limiting the sample to just the smaller subset that is likely to overlap
    # making sure to include the 360/0 degree wrap correctly for WAVES-S!
    if((FieldRAmin>(360-1.6)) | (FieldRAmin<(1.5))){ # at this point we will have the wrap problem to deal with
      Subsel = (FieldRAmin > AngleConversion(Coordinates[1]-1.5)) & (FieldRAmin < AngleConversion(Coordinates[1]+1.5))

    }else{
      # SubSel = ((FieldRAmin-Coordinates[1])<1.5) & ((FieldRAmin-Coordinates[1])>-1.5)
      Subsel = (FieldRAmin > (Coordinates[1]-1.5)) & (FieldRAmin < (Coordinates[1]+1.5))

    }

    if(length(FieldRAmin[SubSel])>0){
      for(ii in 1:length(FieldRAmin[SubSel])){
        if(Coordinates[1]<FieldRAmax[SubSel][ii] & Coordinates[1]>FieldRAmin[SubSel][ii] & Coordinates[2]<FieldDecmax[SubSel][ii] & Coordinates[2]>FieldDecmin[SubSel][ii]){
          # before I accept this, I want to actually check that the object could be on this frame, to account for slight RA/Dec distortions
          # and now determining if the coordinates of this object actually lie within the dimensions of this image
          Coordinates_xy = magWCSradec2xy(c(Coordinates[1]), c(Coordinates[2]), header=ImageHeader[SubSel][[ii]])
          if( (Coordinates_xy[,2]<ImageDimension[SubSel][[ii]]) & (Coordinates_xy[,2]>0) & (Coordinates_xy[,1]<ImageDimension[SubSel][[ii]]) & (Coordinates_xy[,1]>0) ){
            TileNumber[jj]=TileNumber[jj]+1
          }
        }
      }  
    }

  }
  return(TileNumber)
}

magAB2Jansky=function(x){10^(-0.4*(x-8.9))} # taken from ProSpect

duplicateSelection_GAMA=function(Cat, DegreeConfig, ImageDir='/Volumes/ThunderBay/WAVES/profound/tiles_0.3Res/'){
  start_time <- Sys.time()
  cat('beginning Duplicate identification\n')
  
  
  DuplicateTiles=TileNumberCheck(RA_coordinates=Cat$RAmax, Dec_coordinates=Cat$Decmax, Frames=DegreeConfig$frames, 
    FieldRAmin=DegreeConfig$RAminlim, FieldRAmax=DegreeConfig$RAmaxlim, FieldDecmin=DegreeConfig$Decminlim, FieldDecmax=DegreeConfig$Decmaxlim, 
    ImageDir=ImageDir, start_time=start_time)
  CleanSel=DuplicateTiles>1
  Duplicate=rep(1, length(Cat$RAmax))
  Duplicate[!CleanSel]=0
  
  
  print(paste('finished match', (difftime(Sys.time(), start_time, units = "mins")), 'minutes', sep=' '))
  # first, I only want to look in the regions with two overlapping frames
  CleanSel_double=DuplicateTiles==2
  if(length(Cat$RAmax[CleanSel_double])>0){
    # figuring out where the non-finite values are
    print(paste0('duplicate non-finite RA max values: ', length(Cat$RAmax[CleanSel_double][!is.finite(Cat$RAmax[CleanSel_double])]), '/', length(Cat$RAmax[CleanSel_double]) ))
    print(paste0('duplicate non-finite Dec max values: ', length(Cat$Decmax[CleanSel_double][!is.finite(Cat$Decmax[CleanSel_double])]), '/', length(Cat$Decmax[CleanSel_double]) ))
    print(paste0('duplicate non-finite mag values: ', length(Cat$mag[CleanSel_double][!is.finite(Cat$mag[CleanSel_double])]), '/', length(Cat$mag[CleanSel_double]) ))
    clean_max=internalclean(Cat$RAmax[CleanSel_double], Cat$Decmax[CleanSel_double], tiebreak=Cat$mag[CleanSel_double], Nmatch=c(1:10), rad=1) 
                                                                # the nmatch ensures that every object
                                                                # that only appears once is instantly eliminated. 
    print(paste0('duplicate non-finite RA cen values: ', length(Cat$RAcen[CleanSel_double][!is.finite(Cat$RAcen[CleanSel_double])]), '/', length(Cat$RAcen[CleanSel_double]) ))
    print(paste0('duplicate non-finite Dec cen values: ', length(Cat$Deccen[CleanSel_double][!is.finite(Cat$Deccen[CleanSel_double])]), '/', length(Cat$Deccen[CleanSel_double]) ))
    print(paste0('duplicate non-finite mag values: ', length(Cat$mag[CleanSel_double][!is.finite(Cat$mag[CleanSel_double])]), '/', length(Cat$mag[CleanSel_double]) ))
    clean_cen=internalclean(Cat$RAcen[CleanSel_double], Cat$Deccen[CleanSel_double], tiebreak=Cat$mag[CleanSel_double], Nmatch=c(1:10), rad=1) 
                                                              # for some segments, different 
                                                              # versions have resulted in different max coordinate values. 
                                                              # for these scenarios, the cen coordinates will actually be
                                                              # the same, so we want to check fot a duplication of 
                                                              # this coordinate set as well. 
    print('duplicate internal clean completed')                                                          
    Duplicate[CleanSel_double][clean_max]=0
    Duplicate[CleanSel_double][clean_cen]=0
  
  }else{
    print('no duplicates')
  }
  
  # then the regions with three overlapping frames
  CleanSel_triple=DuplicateTiles==3
  if(length(Cat$RAmax[CleanSel_triple])>0){    
    clean_max=internalclean(Cat$RAmax[CleanSel_triple], Cat$Decmax[CleanSel_triple], tiebreak=Cat$mag[CleanSel_triple], Nmatch=c(2:10), rad=1)
    clean_cen=internalclean(Cat$RAcen[CleanSel_triple], Cat$Deccen[CleanSel_triple], tiebreak=Cat$mag[CleanSel_triple], Nmatch=c(2:10), rad=1)
    Duplicate[CleanSel_triple][clean_max]=0
    Duplicate[CleanSel_triple][clean_cen]=0
    print('triplicate internal clean completed') 
  }else{
    print('no triplicates')
  }

  
  # then the regions with four overlapping frames
  CleanSel_quadruple=DuplicateTiles==4
  if(length(Cat$RAmax[CleanSel_quadruple])>0){
    clean_max=internalclean(Cat$RAmax[CleanSel_quadruple], Cat$Decmax[CleanSel_quadruple], tiebreak=Cat$mag[CleanSel_quadruple], Nmatch=c(3:10), rad=1)
    clean_cen=internalclean(Cat$RAcen[CleanSel_quadruple], Cat$Deccen[CleanSel_quadruple], tiebreak=Cat$mag[CleanSel_quadruple], Nmatch=c(3:10), rad=1)
    Duplicate[CleanSel_quadruple][clean_max]=0
    Duplicate[CleanSel_quadruple][clean_cen]=0
    print('quadlicate internal clean completed')      
  }else{
    print('no quadlicates')
  }

return(Duplicate)
}


duplicateSelection_WAVES=function(Cat, DegreeConfig, ImageDir='/Volumes/ThunderBay/WAVES/profound/tiles_0.3Res/'){
  start_time <- Sys.time()
  cat('beginning Duplicate identification\n')
  
  
  # DuplicateTiles=TileNumberCheck(RA_coordinates=Cat$RAmax, Dec_coordinates=Cat$Decmax, Frames=DegreeConfig$frames, 
  #   FieldRAmin=DegreeConfig$RAminlim, FieldRAmax=DegreeConfig$RAmaxlim, FieldDecmin=DegreeConfig$Decminlim, FieldDecmax=DegreeConfig$Decmaxlim, 
  #   ImageDir=ImageDir, start_time=start_time)
  # CleanSel=DuplicateTiles>1
  Duplicate=rep(1, length(Cat$RAmax))
   
  
  # first, I only want to look in the regions with two overlapping frames

  # in principle, I no longer need to worry about the Nmatch criterion for WAVES, as I now want the MOST fragmented solution. 
  clean_max=internalclean(Cat$RAmax, Cat$Decmax, tiebreak=Cat$mag, rad=1, decreasing=TRUE)#, group=TRUE) 
                                                              # decreasing=TRUE ensures that larger values of tiebreak are prioritised, corresponding to smaller flux values. 
                                                              # the impact of this is that the most fragmented solutions will be prioritised over the least fragmented ones. 
  clean_cen=internalclean(Cat$RAcen, Cat$Deccen, tiebreak=Cat$mag, rad=1, decreasing=TRUE)#, group=TRUE) 
                                                            # for some segments, different 
                                                            # versions have resulted in different max coordinate values. 
                                                            # for these scenarios, the cen coordinates will actually be
                                                            # the same, so we want to check for a duplication of 
                                                            # this coordinate set as well. 
  print('duplicate internal clean completed')                                                          
  Duplicate[clean_max]=0
  Duplicate[clean_cen]=0

  print(paste('finished match', (difftime(Sys.time(), start_time, units = "mins")), 'minutes', sep=' '))

return(Duplicate)
}

mask_function = function(catalogue, filelist){
  mask = foreach(file = filelist, .combine='c')%do%{
    bad_poly = fread(file)
    bad_poly[bad_poly$RA > 300,RA := RA - 360]
    return(which(Rwcs_in_poly(catalogue$RAmax, catalogue$Decmax, bad_poly$RA, bad_poly$Dec)))
  }
  
  mask = unique(mask)
  return(mask)
}

duplicate_function = function(catalogue, filelist){
  duplicate = foreach(file = filelist, .combine='c')%do%{
    dup_poly = fread(file)
    dup_poly[dup_poly$RA > 300,RA := RA - 360]
    check = which(Rwcs_in_poly(catalogue$RAmax, catalogue$Decmax, dup_poly$RA, dup_poly$Dec))
    return(check[!catalogue[check,FrameID] %in% names(which.min(table(catalogue[check,FrameID])))])
  }
  
  duplicate = unique(duplicate)
  return(duplicate)
}

# Zpatch = function(Cat){
#   # now computing the approximated Z-band from the i and Y, for regions with missing Z:
#   alpha = 0.4912
#   beta = -0.0281
#   sigma = 0.0134 
#   Zeropoint = 8.9

#   # now adding the equivalent patched version for flux_Zc
#   Cat$mag_Zc = profoundFlux2Mag(Cat$flux_Zc, magzero = Zeropoint)
#   Cat$mag_Yc = profoundFlux2Mag(Cat$flux_Yc, magzero = Zeropoint)
#   Cat$mag_ic = profoundFlux2Mag(Cat$flux_ic, magzero = Zeropoint)
#   Cat$mag_rc = profoundFlux2Mag(Cat$flux_rc, magzero = Zeropoint)

#   Cat$mag_Zc_transform = Cat$mag_Yc + alpha*(Cat$mag_ic - Cat$mag_Yc) + beta 
  
#   # now for the scenario that we have no Y-band data...
#   alpha = -0.7044
#   beta = 0.004
#   sigma = 0.0194 
#   Cat$mag_Zc_transform2 = Cat$mag_ic + alpha*(Cat$mag_rc - Cat$mag_ic) + beta 
  
#   Cat$falseZ = rep(0, length(Cat$uberID))
  
  
#   missingZ = which((Cat$mag_Zt_fake==1) & ((Cat$flux_Zt<0) | is.na(Cat$flux_Zt))) #& ((Cat$mag_Z_transform)<=22))
#   missingZY = which((Cat$mag_Zt_fake==1) & ((Cat$flux_Zt<0) | is.na(Cat$flux_Zt)) & ((Cat$flux_Yt<0) | is.na(Cat$flux_Yt))) #& ((Cat$mag_Z_transform2)<=22)))

#   print(paste0('number of fake Z-band entries:', length(which(Cat$mag_Zt_fake==1))))
#   print(paste0('number of missing Z objects:', length(missingZ)))
#   print(paste0('number of missing Z and Y objects:', length(missingZY)))
  
#   # now actually replacing the relevant values of mag_Zt with the transformed versions, for ease of selection
#   Cat$flux_Zc_patched = Cat$flux_Zc
#   Cat$flux_Zc_patched[missingZ] = profoundMag2Flux(Cat$mag_Zc_transform[missingZ], magzero = Zeropoint) 
#   Cat$flux_Zc_patched[missingZY] = profoundMag2Flux(Cat$mag_Zc_transform2[missingZY], magzero = Zeropoint) 

#   # print(Cat$flux_Zc_patched[missingZ])

#   return(list(flux_Zc_patched = Cat$flux_Zc_patched))
# }

#
#################################################
#
#  Main Code
#

WAVESDir = '/Volumes/ThunderBay/WAVES/'

TilesDir=paste0(WAVESDir, "profound/tiles_0.3Res/")
MeasureDir=paste0(WAVESDir, "profound/measure_0.3Res_WISE/")
PostprocessDir=paste0(WAVESDir, "profound/postprocess_0.3Res_WISE/")
PlotDir=paste0(WAVESDir, "profound/plots/")
RefDir=paste0(WAVESDir, "ref_Sabine/")
CatsDir=paste0(WAVESDir, "cats/")

inputargs=commandArgs(TRUE)
region=as.character(inputargs[1])
version='1p2'
print(paste0('merging catalogues for region ', region))

# all coordinates need to be consistent with https://wavesurvey.org/surveys/
if(region=='WAVES-N'){
  # volume=
  wavesra=c(157.25,157.25,225.0,225.0,157.25)
  wavesdec=c(-3.95,3.95,3.95,-3.95,-3.95)
}else if (region=='WAVES-S'){
  # volume=
  wavesra=c(-30,-30,51.6,51.6,-30) # updated upper RA edge, March 2024
  wavesdec=c(-27,-35.6,-35.6,-27,-27)
}else if (region=='H01'){
  # volume=
  wavesra=c(9.0, 9.0, 21.0, 21.0, 9.0)
  wavesdec=c(-34.0, -29.0, -29.0, -34.0,-34.0)
}else if (region=='H03'){
  # volume=
  wavesra=c(39.0, 39.0, 51.0, 51.0, 39.0)
  wavesdec=c(-34.0, -29.0, -29.0, -34.0,-34.0)
}else if (region=='G09'){
  # volume=
  wavesra=c(129.0, 129.0, 141.0, 141.0, 129.0)
  wavesdec=c(-2.0, 3.0, 3.0, -2.0,-2.0)
}else{
  print('Incorrect region format specified! Needs to be WAVES-N, WAVES-S, H01, H03, G09. ')
}

#
# Read in reference files.
#
gaia=fread(paste0(RefDir,"gaiastarmaskwaves.csv"))
# allr=fread(paste0(RefDir,"allrcounts.csv"))
# tristars=fread(paste0(RefDir,"/ref/",region,"starcounts.csv"))
# x=seq(0,24,0.5)
# xcen=x[1:48]+0.25
# tristarcounts=maghist(tristars$r,breaks=x,plot=FALSE)$counts
DegreeConfig=fread(paste(RefDir,'TileEdges_', region, '.csv', sep=''))
#
# Read in the list of fields to process
#
InputTargetCat=fread(paste0(RefDir, 'Tiles_', region, '.txt'))
#
# For each field read in data
#
if(!file.exists(paste0(CatsDir,region,'_',version,"_initialStitch.parquet"))){
  print('merging catalogue from scratch')
  for (j in 1:length(InputTargetCat$RA)){
    ra=format(round(InputTargetCat$RA[j], 1), nsmall = 1)
    dec=format(round(InputTargetCat$Dec[j], 1), nsmall = 1)
    Filename=paste0(PostprocessDir, "waves_postprocessed_", ra, '_', dec, '.rds')# reading in the postprocessed file
    if(file.exists(Filename)){
      trim=readRDS(Filename)
      datafilex=trim$cat
      cat(date(),Filename,is.data.table(datafilex),"\n")
      #
      if (j==1){datafile0=datafilex}else{datafile0=rbind(datafile0,datafilex)}
    }else{
      print(paste0('WARNING: Missing ', Filename, ' Will need to add output to the initial stitch!!'))
    }
    
  }
  #
  # Determine duplicates from overlap regions
  #
  datafile0$RAcen[!is.finite(datafile0$RAcen)]=datafile0$RAmax[!is.finite(datafile0$RAcen)] # for all NaN cen coordinates, simply adopt the max coordinates
  datafile0$Deccen[!is.finite(datafile0$Deccen)]=datafile0$Decmax[!is.finite(datafile0$Deccen)] # for all NaN cen coordinates, simply adopt the max coordinates
  
  # fwrite(datafile0,file=paste0(CatsDir,region,'_',version,"_initialStitch.csv"),row.names=FALSE) # in case the read-in alone takes a huge amount of time. 
  write_parquet(datafile0, paste0(CatsDir,region,'_',version,"_initialStitch.parquet"))

}else if(!file.exists(paste0(CatsDir,region,'_',version,"_withDuplicate.parquet"))){
  print(paste0('reading in: ', CatsDir,region,'_',version,"_initialStitch.parquet"))
  datafile0=read_parquet(paste0(CatsDir,region,'_',version,"_initialStitch.parquet"))
}

if(!file.exists(paste0(CatsDir,region,'_',version,"_withDuplicate.parquet"))){
  if(version == '1p1' & ((region=='WAVES-S') | (region=='WAVES-N'))){ # now unfortunately a couple of columns went missing in the 1p1 version  
    extraCat = read_parquet(paste0(CatsDir,region,'_1.0_initialStitch.parquet'), 
      col_select=c("uberID", "RAcen","Deccen","RAmax","Decmax","RAGAIA_r","RAGAIA","DecGAIA","DecGAIA_r","RAGAIA_r_cen","DecGAIA_r_cen"))
    datafile0 = datafile0[, (c(119)) := NULL] # removing the extra 'class' column
    print('merging catalogues by uberID')
    datafile0 = merge(datafile0,extraCat,by="uberID")
    
    print('columns in catalogue:')
    print(colnames(datafile0))
  }

  datafile0=datafile0[, 'duplicate' := duplicateSelection_WAVES(datafile0, DegreeConfig, ImageDir=TilesDir)]# Sabine added
  write_parquet(datafile0, paste0(CatsDir,region,'_',version,"_withDuplicate.parquet"))
  print('finished saving parquet')
}else{
  print(paste0('reading in: ', CatsDir,region,'_',version,"_withDuplicate.parquet"))
  datafile0=read_parquet(paste0(CatsDir,region,'_',version,"_withDuplicate.parquet"))
}
# note mag_Zt is already in the catalogue
message('computing mag_rt')
datafile0$mag_rt=8.9-2.5*log10(datafile0$flux_rt)
message('computing mag_it')
datafile0$mag_it=8.9-2.5*log10(datafile0$flux_it)
message('computing mag_Yt')
datafile0$mag_Yt=8.9-2.5*log10(datafile0$flux_Yt)

message('computing mag_rc')
datafile0$mag_rc=8.9-2.5*log10(datafile0$flux_rc)
message('computing mag_ic')
datafile0$mag_ic=8.9-2.5*log10(datafile0$flux_ic)
message('computing mag_Zc')
datafile0$mag_Zc=8.9-2.5*log10(datafile0$flux_Zc)
message('computing mag_Yc')
datafile0$mag_Yc=8.9-2.5*log10(datafile0$flux_Yc)


# datafile0$mag_rt_uncorrected=8.9-2.5*log10(datafile0$flux_rt_uncorrected)

#
# Assign WAVES boundary mask and write out merged catalogue
#
# need to ensure that I'm accounting for the wrapping from 0/30 properly in WAVES-S
if(region=='WAVES-S'){
  datafile0$RAmax[datafile0$RAmax>300] = datafile0$RAmax[datafile0$RAmax>300]-360
}
print('masking WAVES regions')
datafile0[RAmax < min(wavesra) | RAmax > max(wavesra) | Decmax < min(wavesdec) | Decmax > max(wavesdec),"mask"]=1L

print('patching total Z photometry')
#patching the total Z photometry
alpha = 0.4912
beta = -0.0281
sigma = 0.0134

datafile0$mag_Zt_fake = rep(0L, length(datafile0$RAmax))

datafile0=datafile0[is.na(mag_Zt), "mag_Zt_fake" := 1L]
datafile0=datafile0[is.na(mag_Zt), "mag_Zt" := mag_Yt + alpha*(mag_it - mag_Yt) + beta]
datafile0=datafile0[is.na(mag_Zt), "mag_Zc" := mag_Yc + alpha*(mag_ic - mag_Yc) + beta]
datafile0=datafile0[is.na(mag_app_Zt), "mag_app_Zt" := mag_app_Yt + alpha*(mag_app_it - mag_app_Yt) + beta]

alpha = -0.7044
beta = 0.004
sigma = 0.0194

datafile0=datafile0[is.na(mag_Zt), "mag_Zt" := mag_it + alpha*(mag_rt - mag_it) + beta]
datafile0=datafile0[is.na(mag_Zt), "mag_Zc" := mag_ic + alpha*(mag_rc - mag_ic) + beta]
datafile0=datafile0[is.na(mag_app_Zt), "mag_app_Zt" := mag_app_it + alpha*(mag_app_rt - mag_app_it) + beta]

print('masking final galaxies')
# now conducting the final manual masking based on large NGC objects etc
filelist = list.files('/Volumes/ThunderBay/WAVES/NGCmasking_polygons/Masking/', full.names = TRUE)

mask_sel = mask_function(datafile0, filelist)
message(paste0('objects masked: ', length(datafile0$RAmax[mask_sel]) ))
datafile0$mask[mask_sel] = 1L


print('fixing final duplicate flags')
# and now removing the objects that should have been flgged as duplicates in regions where an NGC fix was conducted on an overlapping region
filelist = list.files('/Volumes/ThunderBay/WAVES/NGCmasking_polygons/DuplicateFixing/', full.names = TRUE)


duplicate_sel = duplicate_function(datafile0, filelist)
message(paste0('objects duplicated: ', length(datafile0$RAmax[duplicate_sel]) ))
datafile0$duplicate[duplicate_sel] = 1L


# and then finally removing all sources that are within the region that should be starmasked around a particularly big star that 
# was missed by the GAIA mask (centre (344.415,-29.62), radius 0.55 deg)
if(region=='WAVES-S'){
  print('masking final large star')
  maskedregion=coordmatchsing(RAref=344.415-360,Decref=-29.62,coordcompare=datafile0[,list(RAmax,Decmax)],rad=0.15,inunitref="deg",inunitcompare="deg",radunit="deg")
  maskedobjs=unique(maskedregion$ID[maskedregion$ID>0])
  message(paste0('objects masked for star: ', length(datafile0$RAmax[maskedobjs]) ))
  datafile0$mask[maskedobjs] = 1L

}

print(paste0('number of objects:', length(datafile0$RAmax)))



# removing extra columns
print('removing extra columns')
datafile0 = datafile0[, c('mag_rt', 'mag_it', 'mag_Yt', 'mag_rc', 'mag_ic', 'mag_Yc') := NULL]

print(colnames(datafile0))

# now tranforming back the coordinates, added 30/5/24 in prep for the next version 1p3. 
if(region=='WAVES-S'){
  datafile0$RAmax[datafile0$RAmax< 0] = datafile0$RAmax[datafile0$RAmax < 0]+360
}


print('saving parquet')
write_parquet(datafile0, paste0(CatsDir,region,'_',version,".parquet"))
print('finished saving parquet')
# fwrite(datafile0,file=paste0(CatsDir,region,'_',version,".csv"),row.names=FALSE)
# print('finished saving csv')
# Rfits_write_table(datafile0,file=paste0(CatsDir,region,'_',version,".fits"))
# print('finished saving fits')






# #
# # Generate star-galaxy separation plots
# #
# #   - Generate colour-mag panel
# #
# png(paste0(PlotDir,"stargal_",region,".png"),width=25.0,height=25.0,units="cm",res=240)
# par(mfrow=c(2,1),mar=c(2.0,2.0,2.0,2.0),oma=c(2.0,2.0,0.0,0.0),cex=1.0)
# #
# magplot(datafile0[uberclass=="star",list(mag_rt,mag_Jt-mag_Kt)],pch=".",xlab="r-band magnitude (mag_rt)",ylab="(J-Ks)col",xlim=c(14,24),ylim=c(-1.0,2),col=rgb(0,0,1,0.025))
# points(datafile0[uberclass=="galaxy",list(mag_rt,mag_Jt-mag_Kt)],pch=".",col=rgb(1,0,0,0.025))
# points(datafile0[uberclass=="ambiguous",list(mag_rt,mag_Jt-mag_Kt)],pch=".",col=rgb(0,1,0,0.025))
# legend("topleft",legend=c("Galaxies","Ambiguous","Stars"),col=c(rgb(1,0,0,1),rgb(0,1,0,1),rgb(0,0,1,1)),pch=15,cex=1.5)
# x1=seq(10.0,19.5,0.01)
# y1=0.025+(x1-x1)
# xy1=cbind(x1,y1)
# x2=seq(19.5,30,0.01)
# y2=0.025-0.1*(x2-19.5)^2.0
# xy2=cbind(x2,y2)
# x3=seq(19.5,30,0.01)
# y3=0.025+0.025*(x3-19.5)^1.0
# xy3=cbind(x3,y3)
# xy=rbind(xy1,xy2)
# lines(xy,lty=1,col="black")
# lines(xy3,lty=1,col="black")
# #
# #   - Generate mag v size panel
# #
# magplot(datafile0[uberclass=="star",list(mag_rt,log10(R50))],pch=".",xlab="r-band magnitude (mag_rt)",ylab="R50",xlim=c(14,24),ylim=c(-0.5,0.5),col=rgb(0,0,1,0.025))
# points(datafile0[uberclass=="galaxy",list(mag_rt,log10(R50))],pch=".",col=rgb(1,0,0,0.025))
# points(datafile0[uberclass=="ambiguous",list(mag_rt,log10(R50))],pch=".",col=rgb(0,1,0,0.025))
# legend("topleft",legend=c("Galaxies","Ambiguous","Stars"),col=c(rgb(1,0,0,1),rgb(0,1,0,1),rgb(0,0,1,1)),pch=15,cex=1.5)
# #
# seeing=quantile(datafile0[,log10seeing],prob=0.5)+0.05
# x3=seq(0.0,25.9,0.1)
# y3=seeing-0.075*(x3-20.5)^1.0
# lines(x3,y3,lty=1,col="black",lwd=2)
# x5=seq(20.5,30.0,0.1)
# y5=seeing+(x5-x5)
# lines(x5,y5,lty=1,col="black",lwd=2)
# #
# seeing=quantile(datafile0[,log10seeing],prob=0.05)+0.05
# x3=seq(0.0,25.9,0.1)
# y3=seeing-0.075*(x3-20.5)^1.0
# lines(x3,y3,lty=3,col="black",lwd=2)
# x5=seq(20.5,30.0,0.1)
# y5=seeing+(x5-x5)
# lines(x5,y5,lty=3,col="black",lwd=2)
# #
# seeing=quantile(datafile0[,log10seeing],prob=0.95)+0.05
# x3=seq(0.0,25.9,0.1)
# y3=seeing-0.075*(x3-20.5)^1.0
# lines(x3,y3,lty=3,col="black",lwd=2)
# x5=seq(20.5,30.0,0.1)
# y5=seeing+(x5-x5)
# lines(x5,y5,lty=3,col="black",lwd=2)
# legend("topleft",legend=c("Galaxies","Ambiguous","Stars"),col=c(rgb(1,0,0,1),rgb(0,1,0,1),rgb(0,0,1,1)),pch=15,cex=1.5)
# #
# dev.off()
# #
# # Generate r-band number-count plot
# #
# png(paste0(stub,"/plots/ncounts_",region,".png"),width=15.0,height=15.0,units="cm",res=240)
# #
# magplot(allr[,V1],allr[,V2],log="y",pch=24,col="grey20",cex=0.5,xlab="Total r mag",ylab="Number per 0.5 mag",xlim=c(10,23),ylim=c(0.001,100000))
# magerr(allr[,V1],allr[,V2],yhi=allr[,V3],ylo=allr[,V3])
# lines(xcen,tristarcounts/10.0,col="gold",lwd=3,lty=2)
# #
# x=seq(0,24,0.5)
# xcen=x[1:48]+0.25
# #
# galaxies=maghist(datafile0[(uberclass=="galaxy" & starmask < 1 & mask < 1 & duplicate < 1),mag_rt],breaks=x,plot=FALSE)$counts
# stars=maghist(datafile0[(uberclass=="star" | uberclass=="ambiguous") & starmask < 2 & mask < 1 & duplicate < 1,mag_rt],breaks=x,plot=FALSE)$counts
# cstars=stars
# cstars[1:40]=0.0
# polygon(c(xcen,rev(xcen)),c((galaxies-0.10*cstars-galaxies**0.5)/volume,rev((galaxies+0.10*cstars+galaxies**0.5)/volume)),col=rgb(1,0,0,alpha=0.1),border=rgb(1,0,0,alpha=0.5))
# points(xcen,galaxies/volume,col="red",pch=20,cex=0.5)
# magerr(xcen,galaxies/volume,yhi=(galaxies**0.5)/volume,ylo=(galaxies**0.5)/volume,col="red")
# #
# polygon(c(xcen,rev(xcen)),c(((stars-0.10*cstars+stars**0.5)/volume),rev((stars+0.10*cstars+stars**0.5)/volume)),col=rgb(0,0,1,alpha=0.1),border=rgb(0,0,1,alpha=0.5))
# points(xcen,stars/volume,col="blue",pch=20,cex=0.5)
# magerr(xcen,stars/volume,yhi=(stars**0.5)/volume,ylo=(stars**0.5)/volume,col="blue")
# #
# abline(v=19.8,lty=2,col="black")
# abline(v=21.5,lty=2,col="black")
# #
# ambiguous=maghist(datafile0[class=="ambiguous" & starmask < 2 & mask < 1 & duplicate < 1,mag_rt],breaks=x,plot=FALSE)$counts
# #
# legend("bottomright",legend=c("Galaxies","Stars","Driver et al (2016)","TRILEGAL v1.6"),col=c("red","blue","grey20","gold"),pch=15,cex=1.0)
# #text("topleft",lab=as.character(region),col="black",cex=2.0)
# #
# dev.off()
# #
# png(paste0(stub,"/plots/tile4_",region,".png"),width=20.0,height=10.0,units="cm",res=240)
# par(mar=c(2.0,2.0,2.0,4.0),oma=c(2.0,2.0,0.0,0.0),cex=1.0)
# x=datafile0$RAmax
# y=datafile0$Decmax
# rbPal <- colorRampPalette(c('yellow','orange','darkred','black'))
# z <- rbPal(200)[as.numeric(cut(datafile0$EBV,breaks=200))]
# key <- rbPal(200)[as.numeric(cut(seq(0,0.07,0.07/200.0),breaks=200))]
# magplot(x,y,pch=".",xlim=c(wavesra[1]-0.5,wavesra[4]+0.5),ylim=c(wavesdec[1]-0.5,wavesdec[3]+0.5),col=z,xlab="Right Ascension",ylab="Declination")
# color.legend(wavesra[4]+1.5,wavesdec[1]-0.5,wavesra[4]+1.65,wavesdec[3]+0.5,align="rb",legend=c("0.0","0.01","0.02","0.03","0.04","0.05","0.06","0.07"),rect.col=key,gradient="y")
# mtext(text="E(B-V) (mag)",side=4)
# dev.off()
# #
# datafile0=datafile0[mag_rt < 19.8,]
# #
# # Generate tile layout plot
# #
# png(paste0(stub,"/plots/tile1_",region,".png"),width=20.0,height=10.0,units="cm",res=240)
# par(mar=c(2.0,2.0,2.0,2.0),oma=c(2.0,2.0,0.0,0.0),cex=1.0)
# magplot(datafile0$RAmax,datafile0$Decmax,xlab="Right Ascension",ylab="Declination",xlim=c(wavesra[1]-0.5,wavesra[4]+0.5),ylim=c(wavesdec[1]-0.5,wavesdec[3]+0.5),pch=".",col=hsv(h=120/360,alpha=0.1),side=1:4,labels=c(TRUE,TRUE,FALSE,FALSE))
# dev.off()
# #
# png(paste0(stub,"/plots/tile2_",region,".png"),width=20.0,height=10.0,units="cm",res=240)
# par(mar=c(2.0,2.0,2.0,2.0),oma=c(2.0,2.0,0.0,0.0),cex=1.0)
# magplot(cbind(datafile0[duplicate>0,list(RAmax,Decmax)]),xlab="Right Ascension",ylab="Declination",xlim=c(wavesra[1]-0.5,wavesra[4]+0.5),ylim=c(wavesdec[1]-0.5,wavesdec[3]+0.5),pch=".",col=rgb(0.5,0.5,0.5,alpha=0.1),side=1:4,labels=c(TRUE,TRUE,FALSE,FALSE))
# points(cbind(datafile0[noIR>0,list(RAmax,Decmax)]),pch=".",col=rgb(0,1,1,alpha=0.1))
# points(cbind(datafile0[mask>0,list(RAmax,Decmax)]),pch=".",col=rgb(1,0,0,alpha=0.01))
# points(cbind(datafile0[starmask>0,list(RAmax,Decmax)]),pch=".",col=rgb(0,0,1,alpha=0.1))
# dev.off()
# #
# png(paste0(stub,"/plots/tile3_",region,".png"),width=20.0,height=10.0,units="cm",res=240)
# par(mar=c(2.0,2.0,2.0,2.0),oma=c(2.0,2.0,0.0,0.0),cex=1.0)
# magplot(cbind(datafile0[(class=="galaxy" | class=="ambiguous") & starmask < 1 & mask < 1 & duplicate < 1 & mag_Zt < 21.5,list(RAmax,Decmax)]),xlab="Right Ascension",ylab="Declination",xlim=c(wavesra[1]-0.5,wavesra[4]+0.5),ylim=c(wavesdec[1]-0.5,wavesdec[2]+0.5),pch=".",col=hsv(h=240/360,alpha=0.5),side=1:4,labels=c(TRUE,TRUE,FALSE,FALSE))
# dev.off()
# #
# png(paste0(stub,"/plots/tile5_",region,".png"),width=20.0,height=10.0,units="cm",res=240)
# par(mar=c(2.0,2.0,2.0,2.0),oma=c(2.0,2.0,0.0,0.0),cex=1.0)
# magplot(cbind(datafile0[uberclass=="galaxy" & CATAID==0 & starmask==0 & mask==0 & duplicate==0 & mag_rt < 19.8,list(RAmax,Decmax)]),xlab="Right Ascension",ylab="Declination",xlim=c(wavesra[1]-0.5,wavesra[4]+0.5),ylim=c(wavesdec[1]-0.5,wavesdec[2]+0.5),pch=".",col=hsv(h=300/360,alpha=1.0),side=1:4,labels=c(TRUE,TRUE,FALSE,FALSE))
# length(datafile0[uberclass=="galaxy" & CATAID==0 & starmask==0 & mask==0 & duplicate==0 & mag_rt < 19.8,CATAID])
# dev.off()

