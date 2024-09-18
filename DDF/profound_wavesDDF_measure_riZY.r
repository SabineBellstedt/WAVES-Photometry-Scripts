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
CoreNumber = 2
#
###################
# Setup PATHNAMES #
###################
#

InputDir="/Volumes/ThunderBay/WAVES/profound/tiles_0.3Res_ddf/"
# InputDir2="/Volumes/WAVES/waves/profound/ddfProMosv3/"
DetectDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/detect/'
MeasureDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/measure_WISE/'
PlotDir='/Volumes/ThunderBay/DDF/profound_ddf_v1/plots/'
FixesDir="/Volumes/ThunderBay/DDF/profound_ddf_v1/fix/"

# change filename here to change the list of tiles to be run. 
TileCat=fread('/Volumes/ThunderBay/DDF/profound_ddf_v1/wd03tiles_names.csv')
#
#################################
# reading in the SWarped images #
#################################

registerDoParallel(cores=CoreNumber)

foreach(ii=(1:length(TileCat$ra)))%dopar%{
# foreach(ii=c(13, 9, 6))%dopar%{
# for(ii in c(13, 9, 6)){
  ra = format(round(TileCat$ra[ii], 1), nsmall = 1)
  dec = format(round(TileCat$dec[ii], 1), nsmall = 1)

  message(paste0('*** ', ra, '_', dec, ' ***'))
  
  if(file.exists(paste0(DetectDir,'waves_detect_', ra, '_', dec, '.fits')) & !file.exists(paste0(MeasureDir,'waves_measured_', ra, '_', dec, '.rds'))){

    # Start writing to an output file
    zz <- file(paste0(MeasureDir,'waves_measure_log_', ra, '_', dec, '.txt'), open = "wt")
    sink(zz, type = "output")
    sink(zz, type = "message")
    # sink(file=paste0(MeasureDir,'waves_measure_log_', ra, '_', dec, '.txt'), type='message')
  
    message(paste0('Reading in ProMo images for: ', ra, '_', dec, '\n'))
    
    Available_Images=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
    
    if(file.exists(paste0(InputDir, 'u_', ra, '_', dec, '.fits'))){
      vst_u=Rfits_read_image(paste0(InputDir, 'u_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[1]=TRUE
    }
    if(file.exists(paste0(InputDir, 'g_', ra, '_', dec, '.fits'))){
      vst_g=Rfits_read_image(paste0(InputDir, 'g_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[2]=TRUE
    }
    if(file.exists(paste0(InputDir, 'r_', ra, '_', dec, '.fits'))){
      vst_r=Rfits_read_image(paste0(InputDir, 'r_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[3]=TRUE
    }
    if(file.exists(paste0(InputDir, 'i_', ra, '_', dec, '.fits'))){
      vst_i=Rfits_read_image(paste0(InputDir, 'i_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[4]=TRUE
    }
    if(file.exists(paste0(InputDir, 'Z_', ra, '_', dec, '.fits'))){
      viking_Z=Rfits_read_image(paste0(InputDir, 'Z_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[5]=TRUE
    }
    if(file.exists(paste0(InputDir, 'Y_', ra, '_', dec, '.fits'))){
      viking_Y=Rfits_read_image(paste0(InputDir, 'Y_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[6]=TRUE
    }
    if(file.exists(paste0(InputDir, 'J_', ra, '_', dec, '.fits'))){
      viking_J=Rfits_read_image(paste0(InputDir, 'J_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[7]=TRUE
    }
    if(file.exists(paste0(InputDir, 'H_', ra, '_', dec, '.fits'))){
      viking_H=Rfits_read_image(paste0(InputDir, 'H_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[8]=TRUE
    }
    if(file.exists(paste0(InputDir, 'Ks_', ra, '_', dec, '.fits'))){
      viking_Ks=Rfits_read_image(paste0(InputDir, 'Ks_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[9]=TRUE
    }
    if(file.exists(paste0(InputDir, 'w1_', ra, '_', dec, '.fits'))){
      wise_W1=Rfits_read_image(paste0(InputDir, 'w1_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[10]=TRUE
    }
    if(file.exists(paste0(InputDir, 'w2_', ra, '_', dec, '.fits'))){
      wise_W2=Rfits_read_image(paste0(InputDir, 'w2_', ra, '_', dec, '.fits'), ext=1)
      Available_Images[11]=TRUE
    }
    
    # for those images that are missing, replicate the image in an existing band
    index=1:9
    FirstAvailableImage=index[Available_Images>0][1]
    if(FirstAvailableImage==1){
      DummyImage = vst_u
    }else if(FirstAvailableImage==2){
      DummyImage = vst_g
    }else if(FirstAvailableImage==3){
      DummyImage = vst_r
    }else if(FirstAvailableImage==4){
      DummyImage = vst_i
    }else if(FirstAvailableImage==5){
      DummyImage = viking_Z
    }else if(FirstAvailableImage==6){
      DummyImage = viking_Y
    }else if(FirstAvailableImage==7){
      DummyImage = viking_J
    }else if(FirstAvailableImage==8){
      DummyImage = viking_H
    }else if(FirstAvailableImage==9){
      DummyImage = viking_Ks
    }else if(FirstAvailableImage==10){
      DummyImage = wise_W1
    }else if(FirstAvailableImage==11){
      DummyImage = wise_W2
    }
    
    if(Available_Images[1]==FALSE){
      vst_u = DummyImage
      vst_u$imDat[]=NA
      message('*** Missing image in the u band. Generating all-NA image ***')
    }
    if(Available_Images[2]==FALSE){
      vst_g = DummyImage
      vst_g$imDat[]=NA
      message('*** Missing image in the g band. Generating all-NA image ***')
    }
    if(Available_Images[3]==FALSE){
      vst_r = DummyImage
      vst_r$imDat[]=NA
      message('*** Missing image in the r band. Generating all-NA image ***')
    }
    if(Available_Images[4]==FALSE){
      vst_i = DummyImage
      vst_i$imDat[]=NA
      message('*** Missing image in the i band. Generating all-NA image ***')
    }
    if(Available_Images[5]==FALSE){
      viking_Z = DummyImage
      viking_Z$imDat[]=NA
      message('*** Missing image in the Z band. Generating all-NA image ***')
    }
    if(Available_Images[6]==FALSE){
      viking_Y = DummyImage
      viking_Y$imDat[]=NA
      message('*** Missing image in the Y band. Generating all-NA image ***')
    }
    if(Available_Images[7]==FALSE){
      viking_J = DummyImage
      viking_J$imDat[]=NA
      message('*** Missing image in the J band. Generating all-NA image ***')
    }
    if(Available_Images[8]==FALSE){
      viking_H = DummyImage
      viking_H$imDat[]=NA
      message('*** Missing image in the H band. Generating all-NA image ***')
    }
    if(Available_Images[9]==FALSE){
      viking_Ks = DummyImage
      viking_Ks$imDat[]=NA
      message('*** Missing image in the Ks band. Generating all-NA image ***')
    }
    if(Available_Images[10]==FALSE){
      wise_W1 = DummyImage
      wise_W1$imDat[]=NA
      message('*** Missing image in the W1 band. Generating all-NA image ***')
    }
    if(Available_Images[11]==FALSE){
      wise_W2 = DummyImage
      wise_W2$imDat[]=NA
      message('*** Missing image in the W2 band. Generating all-NA image ***')
    }
    
    #######################
    # Run ProFound Script #
    #######################
    #
    detect_image=Rfits_read_image(paste0(DetectDir,'waves_detect_', ra, '_', dec, '.fits'), header=FALSE, ext='image')
    detect_segim=Rfits_read_image(paste0(DetectDir,'waves_detect_', ra, '_', dec, '.fits'), header=FALSE, ext='segim_orig')
    dilated_segim=Rfits_read_image(paste0(DetectDir,'waves_detect_', ra, '_', dec, '.fits'), header=FALSE, ext='segim')
    # detect_skyRMS=Rfits_read_image(paste0(DetectDir,'waves_detect_', ra, '_', dec, '.fits'), header=FALSE, ext='skyRMS')
    # use the segim_orig instead of segim_dilate so that we can output undilated colour photometry. 
    if(file.exists(paste0(FixesDir,'waves_segID_merge_', ra, '_', dec, '.rds'))){
      segID_merge_new=readRDS(paste0(FixesDir,'waves_segID_merge_', ra, '_', dec, '.rds'))
      fixed_segim=profoundSegimKeep(detect_segim, segID_merge = segID_merge_new)
    }else if(file.exists(paste0(FixesDir,'waves_segID_merge_auto_', ra, '_', dec, '.rds'))){
      segID_merge_new=readRDS(paste0(FixesDir,'waves_segID_merge_auto_', ra, '_', dec, '.rds'))$segID_merge
      fixed_segim=profoundSegimKeep(detect_segim, segID_merge = segID_merge_new)
      message(paste0('**** Using automatic fixing for ', ra, '_', dec, '. Manual fixing not completed ****'))
    }else{ # in this scenario we don't bother fixing at all

      message(paste0('**** No fixing available for ', ra, '_', dec, '. ****'))
      fixed_segim=detect_segim
    }
    
    #
    InputImages=list(vst_u, vst_g, vst_r, vst_i, viking_Z, viking_Y, viking_J, viking_H, viking_Ks, wise_W1, wise_W2)
    
    everything=profoundMultiBand(
      inputlist = InputImages,
      segim = fixed_segim, 
      # skyRMS = detect_skyRMS, 
      skycut = 2.0, 
      pixcut = 13,
      ext = 1,
      tolerance = 15,
      reltol = -10,
      cliptol = 100,
      iters_det = 6,
      iters_tot = c(2,2,2,2,2,2,2,2,2,3,3), 
      totappend='t',
      sizes_tot = c(5,5,5,5,5,5,5,5,5,15,15),
      colappend = 'c',
      detectbands = 'Z',
      multibands = c('u','g','r','i','Z','Y','J','H','Ks','W1','W2'),
      keepsegims = TRUE,
      magzero = c(23.9,23.9,23.9,23.9,23.9,23.9,23.9,23.9,23.9,23.183,22.819), # v2 images are in microJy (not sure about WISE though!)
      dotot = TRUE,
      docol = TRUE,
      dogrp = TRUE,
      verbose = TRUE,
      box = c(100,100,100,100,100,100,100,100,100,200,200),
      boxiters = 4,
      boxadd = 50,
      grid = 50,
      roughpedestal = TRUE,
      redosegim = FALSE,
      deblend = FALSE,
      groupstats = TRUE,
      mask = 0,
      SBdilate = 1.0,
      SBN100 = 100,
      app_diam = 1.4,#4MOST fibre aperture
      fluxtype = 'Jansky'
      # redosky=FALSE
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
    CairoPDF(file=paste0(PlotDir,'diagnostics_measure_', ra, '_', dec, '.pdf'),width=24.0,height=24.0)
    plot(everything$pro_detect)
    dev.off()
    ###

    # Stop writing to the output file
    sink()
  
  }  
} 