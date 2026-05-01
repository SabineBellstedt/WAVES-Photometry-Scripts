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
library(yaml)

inputargs=commandArgs(TRUE)
configFilename=as.character(inputargs[1])
config = yaml.load_file(configFilename)

MeasureDir = paste0(config$path$measure, config$version$detect, config$version$measure, '/')
PostprocessDir = paste0(config$path$postprocess, config$version$detect, config$version$measure, config$version$postprocess, '/')
CatsDir = paste0(config$path$catalogues, config$version$detect, config$version$measure, config$version$postprocess, config$version$catalogue, '/')
PlotDir = paste0(CatsDir, "plots/")

MeasureVersionSuffix = paste0('_', config$version$detect, config$version$measure)
PostprocessVersionSuffix = paste0('_', config$version$detect, config$version$measure, config$version$postprocess)
CatalogueVersionSuffix = paste0('_', config$version$detect, config$version$measure, config$version$postprocess, config$version$catalogue)

print(CatalogueVersionSuffix)

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

magAB2Jansky=function(x){10^(-0.4*(x-8.9))} # taken from ProSpect


duplicateSelection_WAVES=function(Cat){
  start_time <- Sys.time()
  cat('beginning Duplicate identification\n')
  
  Duplicate=rep(1, length(Cat$RAmax))

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

#
#################################################
#
#  Main Code
# 
# Read in reference files.
#
gaia=fread(paste0(config$path$reference, config$referencefiles$gaiastars))

for(region in config$regions){
  InputTargetCat=fread(paste0(config$path$reference, config$referencefiles$tilelist[[region]]))
  #
  # For each field read in data
  #
  if(!file.exists(paste0(CatsDir,region,CatalogueVersionSuffix,"_initialStitch.parquet"))){
    print('merging catalogue from scratch')
    for (j in 1:length(InputTargetCat$RA)){
      ra=format(round(InputTargetCat$RA[j], 1), nsmall = 1)
      dec=format(round(InputTargetCat$Dec[j], 1), nsmall = 1)
      PostprocessFilename=paste0(PostprocessDir, "waves_postprocessed_", ra, '_', dec, PostprocessVersionSuffix, '.parquet')# reading in the postprocessed file
      if(file.exists(PostprocessFilename)){
        # trim=readRDS(PostprocessFilename)
        # datafilex=trim$cat
        datafilex = read_parquet(PostprocessFilename)
        cat(date(),PostprocessFilename,is.data.table(datafilex),"\n")
        #
        if (j==1){datafile0=datafilex}else{datafile0=rbind(datafile0,datafilex)}
      }else{
        print(paste0('WARNING: Missing ', PostprocessFilename, ' Will need to add output to the initial stitch!!'))
      }     
    }
    #
    # Determine duplicates from overlap regions
    #
    datafile0$RAcen[!is.finite(datafile0$RAcen)]=datafile0$RAmax[!is.finite(datafile0$RAcen)] # for all NaN cen coordinates, simply adopt the max coordinates
    datafile0$Deccen[!is.finite(datafile0$Deccen)]=datafile0$Decmax[!is.finite(datafile0$Deccen)] # for all NaN cen coordinates, simply adopt the max coordinates
    
    write_parquet(datafile0, paste0(CatsDir,region,CatalogueVersionSuffix,"_initialStitch.parquet"))
  
  }else if(!file.exists(paste0(CatsDir,region,CatalogueVersionSuffix,"_withDuplicate.parquet"))){
    print(paste0('reading in: ', CatsDir,region,CatalogueVersionSuffix,"_initialStitch.parquet"))
    datafile0=read_parquet(paste0(CatsDir,region,CatalogueVersionSuffix,"_initialStitch.parquet"))
  }
  
  if(!file.exists(paste0(CatsDir,region,CatalogueVersionSuffix,"_withDuplicate.parquet"))){  
    datafile0=datafile0[, 'duplicate' := duplicateSelection_WAVES(datafile0)]# Sabine added
    # write_parquet(datafile0, paste0(CatsDir,region,CatalogueVersionSuffix,"_withDuplicate.parquet"))
    # print('finished saving parquet')
  }else{
    print(paste0('reading in: ', CatsDir,region,CatalogueVersionSuffix,"_withDuplicate.parquet"))
    datafile0=read_parquet(paste0(CatsDir,region,CatalogueVersionSuffix,"_withDuplicate.parquet"))
  }
  
  #
  # Assign WAVES boundary mask and write out merged catalogue
  #
  # need to ensure that I'm accounting for the wrapping from 0/30 properly in WAVES-S
  if(region=='WAVES-S'){
    datafile0$RAmax[datafile0$RAmax>300] = datafile0$RAmax[datafile0$RAmax>300]-360
  }
  print('masking WAVES regions')
  # datafile0[RAmax < min(wavesra) | RAmax > max(wavesra) | Decmax < min(wavesdec) | Decmax > max(wavesdec),"mask"]=1L
  datafile0[RAmax < config$regioncoordinates[[region]]$RA[1] | 
    RAmax > config$regioncoordinates[[region]]$RA[2] | 
    Decmax < config$regioncoordinates[[region]]$Dec[1] | 
    Decmax > config$regioncoordinates[[region]]$Dec[2],"mask"]=1L 

  print('fixing final duplicate flags')
  # and now removing the objects that should have been flgged as duplicates in regions where an NGC fix was conducted on an overlapping region
  filelist = list.files(config$referencefiles$duplicatefixpolygons, full.names = TRUE)
    
  duplicate_sel = duplicate_function(datafile0, filelist)
  message(paste0('objects duplicated: ', length(datafile0$RAmax[duplicate_sel]) ))
  datafile0$duplicate[duplicate_sel] = 1L
  
  
  print(paste0('number of objects:', length(datafile0$RAmax)))
  
  # removing extra columns (these can all be regenerated from flux columns)
  print('removing extra columns')
  datafile0 = datafile0[, c('mag_rt', 'mag_it', 'mag_Yt', 'mag_rc', 'mag_ic', 'mag_Yc') := NULL]
  
  if(region=='WAVES-S'){
    datafile0$RAmax[datafile0$RAmax< 0] = datafile0$RAmax[datafile0$RAmax < 0]+360
  }
  
  
  print('saving parquet')
  write_parquet(datafile0, paste0(CatsDir,region,CatalogueVersionSuffix,".parquet"))
  print('finished saving parquet')

}