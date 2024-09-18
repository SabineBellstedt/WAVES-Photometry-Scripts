# .libPaths(c('/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library/',.libPaths()))
library(celestial)
# library(devtools)
#
library(Cairo)
library(ProFound) 
library(magicaxis)
library(data.table)
library(Rcpp)
require(foreign)
require(MASS)
library(Rfits)
library(doParallel)
library(bit64)
CoreNumber = 4
#
###################
# Setup PATHNAMES #
###################
#
# inputargs=commandArgs(TRUE)
# ra=as.character(inputargs[1])
# dec=as.character(inputargs[2])

region = 'G09'


WAVESwideDir = '/Volumes/ThunderBay/WAVES/'

InputDir = paste0(WAVESwideDir, "profound/tiles_0.3Res/")
DetectDir = paste0(WAVESwideDir, "profound/detect_0.3Res/")
MeasureDir = paste0(WAVESwideDir, "profound/measure_0.3Res_WISE/")
PlotDir = paste0(WAVESwideDir, "profound/plots/")
FixesDir = paste0(WAVESwideDir, "profound/segmentFixes/")
RefDir=paste0(WAVESwideDir, "ref_Sabine/")

InputTargetCat=fread(paste0(RefDir, 'Tiles_', region, '.txt'))

#
#################################
# reading in the SWarped images #
#################################

# only bother if the output doesn't already exist
registerDoParallel(cores=CoreNumber)

foreach(j=(1:length(InputTargetCat$RA)))%dopar%{

  ra=as.character((format(round(InputTargetCat$RA[j], 1), nsmall = 1)))
  dec=as.character((format(round(InputTargetCat$Dec[j], 1), nsmall = 1)))


  if(!file.exists(paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds'))){
    print(paste0('Reading in SWarped images for: ', ra, '_', dec))
    
    #
    kids_u=Rfits_point(paste0(InputDir, 'u_', ra, '_', dec, '.fits'), ext=1)
    kids_g=Rfits_point(paste0(InputDir, 'g_', ra, '_', dec, '.fits'), ext=1)
    kids_r=Rfits_point(paste0(InputDir, 'r_', ra, '_', dec, '.fits'), ext=1)
    kids_i=Rfits_point(paste0(InputDir, 'i_', ra, '_', dec, '.fits'), ext=1)
    viking_Z=Rfits_point(paste0(InputDir, 'Z_', ra, '_', dec, '.fits'), ext=1)
    viking_Y=Rfits_point(paste0(InputDir, 'Y_', ra, '_', dec, '.fits'), ext=1)
    viking_J=Rfits_point(paste0(InputDir, 'J_', ra, '_', dec, '.fits'), ext=1)
    viking_H=Rfits_point(paste0(InputDir, 'H_', ra, '_', dec, '.fits'), ext=1)
    viking_Ks=Rfits_point(paste0(InputDir, 'Ks_', ra, '_', dec, '.fits'), ext=1)
    WISE_1=Rfits_point(paste0(InputDir, 'w1_', ra, '_', dec, '.fits'), ext=1)
    WISE_2=Rfits_point(paste0(InputDir, 'w2_', ra, '_', dec, '.fits'), ext=1)
    
    
    #######################
    # Run ProFound Script #
    #######################
    #
    detect_segim=Rfits_read(paste0(DetectDir,'waves_detect_', ra, '_', dec, '.fits'), pointer=FALSE, header=FALSE)$segim_orig
    dilated_segim=Rfits_read(paste0(DetectDir,'waves_detect_', ra, '_', dec, '.fits'), pointer=FALSE, header=FALSE)$segim
    segID_merge_new=readRDS(paste0(FixesDir,'waves_segID_merge_', ra, '_', dec, '.rds'))
    # use the segim_orig instead of segim_dilate so that we can output undilated colour photometry. 
    fixed_segim=profoundSegimKeep(detect_segim, segID_merge = segID_merge_new)
    #
    InputImages=list(kids_u, kids_g, kids_r, kids_i, viking_Z, viking_Y, viking_J, viking_H, viking_Ks, WISE_1, WISE_2)
    
    everything=profoundMultiBand(
      inputlist=InputImages, # proMB
      segim=fixed_segim,  # proMB
      skycut = 2.0,  # pro
      pixcut = 13, # pro
      ext = 1, # pro
      tolerance = 15, # pro
      reltol = -10, # pro
      cliptol = 100, # pro
      iters_det = 6, # proMB
      iters_tot = c(2,2,2,2,2,2,2,2,2,3,3),  # proMB
      totappend='t', # proMB
      sizes_tot = c(5,5,5,5,5,5,5,5,5,15,15), # proMB
      colappend='c', # proMB
      detectbands=c('r','i','Z','Y'), # proMB
      multibands=c('u','g','r','i','Z','Y','J','H','K','W1','W2'), # proMB
      keepsegims = TRUE, # proMB
      magzero = c(0,0,0,0,30,30,30,30,30,23.183,22.819), # proMB
      dotot = TRUE, # proMB
      docol = TRUE, # proMB
      dogrp = TRUE, # proMB
      verbose = TRUE, # proMB
      box = c(100,100,100,100,100,100,100,100,100,200,200), # pro
      boxiters = 4, # pro
      boxadd = c(50,50,50,50,50,50,50,50,50,50,50), # pro
      grid = c(50,50,50,50,50,50,50,50,50,50,50), # pro
      roughpedestal = TRUE, # pro
      redosegim = FALSE, # pro
      deblend = FALSE, # pro
      groupstats = TRUE, # pro
      mask = 0, # pro
      SBdilate = 1.0, # pro
      SBN100 = 100, # pro
      app_diam = 1.4,#4MOST fibre aperture,  # pro
      fluxtype = 'Jansky'  # pro
    )
    
    
    #
    # save file
    #
    everything$pro_detect$fixed_segim=fixed_segim
    everything$pro_detect$detect_segim=detect_segim
    everything$pro_detect$dilated_segim=dilated_segim
    everything$pro_detect$segID_merge=segID_merge_new
    saveRDS(everything,file=paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds'))
    #
    CairoPDF(file=paste0(PlotDir,'diagnostics_wise_', ra, '_', dec, '.pdf'),width=24.0,height=24.0)
    plot(everything$pro_detect)
    dev.off()
    ###
  }else{
    print(paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds already exists'))
  }
}