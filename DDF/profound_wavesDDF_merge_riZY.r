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
library(sphereplot)
# library(fst)
library(arrow)
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


duplicateSelection_WAVES=function(Cat){
  start_time <- Sys.time()
  cat('beginning Duplicate identification\n')
  cat(paste0('length of catalogue: ', length(Cat$RAmax), '\n'))
  
  
  Duplicate=rep(1, length(Cat$RAmax))
   
  
  # first, I only want to look in the regions with two overlapping frames

  # in principle, I no longer need to worry about the Nmatch criterion for WAVES, as I now want the MOST fragmented solution. 
  cat('max duplicates\n')
  clean_max=internalclean(Cat$RAmax, Cat$Decmax, tiebreak=Cat$mag, rad=1, decreasing=TRUE, iter=TRUE) 
                                                              # decreasing=TRUE ensures that larger values of tiebreak are prioritised, corresponding to smaller flux values. 
                                                              # the impact of this is that the most fragmented solutions will be prioritised over the least fragmented ones. 
                                                              # turning on iter=TRUE to ensure no remaining targets are within the separation.                                                      
  # cat('cen duplicates\n')
  # clean_cen=internalclean(Cat$RAcen, Cat$Deccen, tiebreak=Cat$mag, rad=1, decreasing=TRUE, iter=TRUE) 
                                                            # for some segments, different 
                                                            # versions have resulted in different max coordinate values. 
                                                            # for these scenarios, the cen coordinates will actually be
                                                            # the same, so we want to check for a duplication of 
                                                            # this coordinate set as well. 
  print('duplicate internal clean completed')                                                          
  Duplicate[clean_max]=0
  # Duplicate[clean_cen]=0

  print(paste('finished match', (difftime(Sys.time(), start_time, units = "mins")), 'minutes', sep=' '))

return(Duplicate)
}

#
#################################################
#
#  Main Code
#

RefDir="/Volumes/ThunderBay/DDF/profound_ddf_v0/ref_Sabine/"


MeasureDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/measure_WISE/'
PostprocessDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/postprocess_WISE_zpTest_GAAPoffset/'
PlotDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/plots/'
RefDir="/Volumes/ThunderBay/DDF/profound_ddf_v0/ref_Sabine/"
CatsDir="/Volumes/ThunderBay/DDF/profound_ddf_v1/"

# region='WD10'
version='2p4_zpTest4' # the third version
# for(region in c('WD01', 'WD02', 'WD03', 'WD10')){
for(region in c('WD01')){

  print(paste0('merging catalogues for region ', region))
  
  
  if(region=='WD01'){
    # wavesra = 9.475 # 00:37:54
    # wavesdec = -43.95 # -43:57:00
    wavesra = 9.5 # 00:38:00
    wavesdec =  -43.95 # -43:57:00
    aesop_rotation = 0
  }else if(region=='WD02'){
    # wavesra = 35.708
    # wavesra = 35.865 # possible shift to move out of starmask
    # wavesdec = -5.025
    wavesra = 35.875 # 02:23:30
    wavesdec = -5.025 # -5:01:30
    aesop_rotation = 30
  }else if(region=='WD03'){
    # wavesra = 53.125
    # wavesdec = -28.1
    wavesra = 53.125 # 03:32:30
    wavesdec = -28.1 # -28:06:00
    aesop_rotation = 30
  }else if(region=='WD10'){
    # wavesra = 150.1
    # wavesdec = 2.182
    wavesra = 150.125 # 10:00:30
    wavesdec = 2.2 # +02:12:00
    aesop_rotation = 30
  }else{
    print('Incorrect region format specified! Needs to be WD01, WD02, WD03, WD10. ')
  }
  
  
  #
  # Read in reference files.
  #
  gaia=fread(paste0(RefDir,"gaiastarmaskwaves.csv"))
  
  #
  # Read in the list of fields to process
  #
  InputTargetCat=paste0(RefDir, 'ddftiles_names.csv')
  
  Cat=fread(InputTargetCat)
  
  Buffer=0.5
  
  if(region=='WD01'){
    RegionSel = (Cat$ra>(6-Buffer)) & (Cat$ra<(12+Buffer)) & (Cat$dec>(-46-Buffer)) & (Cat$dec< (-42+Buffer))
  }else if(region=='WD02'){
    RegionSel = (Cat$ra>(33-Buffer)) & (Cat$ra<(38+Buffer)) & (Cat$dec>(-7-Buffer)) & (Cat$dec< (-3+Buffer))
  }else if(region=='WD03'){
    RegionSel = (Cat$ra>(51-Buffer)) & (Cat$ra<(56+Buffer)) & (Cat$dec>(-30-Buffer)) & (Cat$dec< (-26+Buffer))
  }else if(region=='WD10'){
    RegionSel = (Cat$ra>(148-Buffer)) & (Cat$ra<(152+Buffer)) & (Cat$dec>(0-Buffer)) & (Cat$dec< (4+Buffer))
  }
  
  InputTargetCat=Cat[RegionSel]
  #
  # For each field read in data
  #
  if(file.exists(paste0(CatsDir,region,'_init_',version,".parquet"))){
    datafile0 = read_parquet(paste0(CatsDir,region,'_init_',version,".parquet"))
  }else{
    FirstFile=TRUE
    print('merging catalogue from scratch')
    for (j in 1:length(InputTargetCat$ra)){
      ra=format(round(InputTargetCat$ra[j], 1), nsmall = 1)
      dec=format(round(InputTargetCat$dec[j], 1), nsmall = 1)
      Filename=paste0(PostprocessDir, "waves_postprocessed_", ra, '_', dec, '.rds')# reading in the postprocessed file
      if(file.exists(Filename)){
        print(Filename)
        trim=readRDS(Filename)
        datafilex=trim$cat
        cat(date(),Filename,is.data.table(datafilex),"\n")
        #
        if (FirstFile){datafile0=datafilex}else{datafile0=rbind(datafile0,datafilex)} 
        FirstFile=FALSE   
      }else{
        message(paste0(Filename, ' does not exist!'))
      }
    
    }
    #
    # Determine duplicates from overlap regions
    #
    datafile0$RAcen[!is.finite(datafile0$RAcen)]=datafile0$RAmax[!is.finite(datafile0$RAcen)] # for all NaN cen coordinates, simply adopt the max coordinates
    datafile0$Deccen[!is.finite(datafile0$Deccen)]=datafile0$Decmax[!is.finite(datafile0$Deccen)] # for all NaN cen coordinates, simply adopt the max coordinates
  
    write_parquet(datafile0,paste0(CatsDir,region,'_init_',version,".parquet"))
    print('finished saving initial parquet')
  
  }

  # checking for duplicated uberIDs
  print('duplicated uberIDs?')
  print(anyDuplicated(datafile0$uberID))
  
  
  datafile0=datafile0[, 'duplicate' := duplicateSelection_WAVES(datafile0)]
  #
  datafile0$mag_rt=8.9-2.5*log10(datafile0$flux_rt)
  datafile0$mag_it=8.9-2.5*log10(datafile0$flux_it)
  datafile0$mag_Yt=8.9-2.5*log10(datafile0$flux_Yt)
  
  datafile0$mag_rc=8.9-2.5*log10(datafile0$flux_rc)
  datafile0$mag_ic=8.9-2.5*log10(datafile0$flux_ic)
  datafile0$mag_Zc=8.9-2.5*log10(datafile0$flux_Zc)
  datafile0$mag_Yc=8.9-2.5*log10(datafile0$flux_Yc)
  
  
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
  
  #
  # Assign WAVES boundary mask and write out merged catalogue
  #
  #
  # Define and rotate 4MOST footprints
  #
  p2s = 60e-3 #pos in mm -> coord in asec
  AESOP_layout = Rfits_read(paste0(RefDir, "fs_layout_20210824.fits"))
  AESOP_field_deg = AESOP_layout[[2]][,2:3]/p2s/3600
  #
  temp = sph2car(AESOP_field_deg)
  temp = rotate3d(temp, x=1, y=0, z=0, angle=aesop_rotation*pi/180)
  aesop=car2sph(temp)
  aesopx=aesop[,1]
  aesopx[31]=aesopx[1]
  aesopy=aesop[,2]
  aesopy[31]=aesopy[1]
  #
  AESOP_tight_hex = cbind(x = c(0, 1.097046, 1.097046, 1.097046, 0, -1.097046, -1.097046, -1.097046, 0),
                          y = c(1.26676, 0.6333799, 0, -0.6333799, -1.26676, -0.6333799, 0, 0.6333799, 1.26676)
  )
  #
  AESOP_tight_circle = cbind(x = sin(0:359*pi/180) * 1.097046,
                             y = cos(0:359*pi/180) * 1.097046
  )
  temp = sph2car(AESOP_tight_hex)
  temp = rotate3d(temp, x=1, y=0, z=0, angle=aesop_rotation*pi/180)
  aesopt=car2sph(temp)
  aesoptx=aesopt[,1]
  aesoptx[10]=aesoptx[1]
  aesopty=aesopt[,2]
  aesopty[10]=aesopty[1]
  
  ddfra=wavesra+aesoptx[c(1,2,4,5,6,8,9)]/cosd(wavesdec+aesopty[c(1,2,4,5,6,8,9)])
  ddfdec=wavesdec+aesopty[c(1,2,4,5,6,8,9)]
  
  inhex=profoundInPoly(x=datafile0$RAmax,y=datafile0$Decmax,poly_x=ddfra,poly_y=ddfdec)
  datafile0[!inhex,"mask"]=1L

  print(ddfra)
  print(ddfdec)


  # # finally resaving the columns with Ks to be K
  # colnames(datafile0)[colnames(datafile0) == "mag_app_Kst"] ="mag_app_Kt"
  # colnames(datafile0)[colnames(datafile0) == "mag_Kst"] ="mag_Kt"
  # colnames(datafile0)[colnames(datafile0) == "flux_Kst"] ="flux_Kt"
  # colnames(datafile0)[colnames(datafile0) == "flux_err_Kst"] ="flux_err_Kt"
  # colnames(datafile0)[colnames(datafile0) == "flux_Kst_uncorrected"] ="flux_Kt_uncorrected"
  # colnames(datafile0)[colnames(datafile0) == "flux_err_Kst_uncorrected"] ="flux_err_Kt_uncorrected"
  # colnames(datafile0)[colnames(datafile0) == "flux_Ksc"] ="flux_Kc"
  # colnames(datafile0)[colnames(datafile0) == "flux_err_Ksc"] ="flux_err_Kc"
  # colnames(datafile0)[colnames(datafile0) == "flux_Ksc_uncorrected"] ="flux_Kc_uncorrected"
  # colnames(datafile0)[colnames(datafile0) == "flux_err_Ksc_uncorrected"] ="flux_err_Kc_uncorrected"

  # removing extra columns
  print('removing extra columns')
  datafile0 = datafile0[, c('mag_rt', 'mag_it', 'mag_Yt', 'mag_rc', 'mag_ic', 'mag_Yc') := NULL]

  datafile0$uberID = as.character(datafile0$uberID)
  
  print(colnames(datafile0))

  
  write_parquet(datafile0, paste0(CatsDir,region,'_',version,".parquet"))
  print('finished saving parquet')
  

  

}  
  