# .libPaths(c('/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library/',.libPaths()))
library(data.table)
library(ProFound)
library(magicaxis)
library(Rfits)
library(celestial)
library(Rwcs)
library(doParallel)
CoreNumber = 1


TilesDir = "/Volumes/WAVES/waves/profound/ddfProMosv2/" # specify location of saved detection segment files
TilesDir2 = "/Volumes/WAVES/waves/profound/ddfProMosv3/" # specify location of saved detection segment files
DataDir = '/Volumes/ThunderBay/DDF/profound_ddf_v1/detect/' # specify location of saved detection segment files
OutputDir="/Volumes/ThunderBay/DDF/profound_ddf_v1/fix/"
gaia=fread("/Volumes/WAVES/waves/ref/GAIA/gaiastarmaskwaves.csv")


TileCat=fread('/Volumes/ThunderBay/DDF/profound_ddf_v1/wd10tiles_names.csv')
TileCat$Frame = paste0(trimws(format(round(TileCat$ra, 1), nsmall = 1)), '_', trimws(format(round(TileCat$dec, 1), nsmall = 1)))

# in the first run of the DDF catalogues, we won't have the KiDS data yet, so I will want to try this script using just Z and Y somehow...
registerDoParallel(cores=CoreNumber)

foreach(frameID=1:length(TileCat$Frame))%dopar%{
# for(frameID in 1:length(TileCat$Frame)){
  frame=TileCat$Frame[frameID]
  segimFilename_Input=paste(DataDir, 'waves_detect_', frame, '.fits', sep='')
  if(file.exists(segimFilename_Input) & !file.exists(paste0(OutputDir, 'waves_segID_merge_auto_', frame, '.rds'))){
    cat(paste0('reading in frame ', frame, ' \n') )
    SegimCat=Rfits_read(segimFilename_Input, pointer=FALSE, header=FALSE)
    cat('fits file read in \n')
  
    gimage=Rfits_read_image(paste0(TilesDir, 'g_', frame,'.fits'), ext=1)
    rimage=Rfits_read_image(paste0(TilesDir, 'r_', frame,'.fits'), ext=1)     
    Zimage=Rfits_read_image(paste0(TilesDir2, 'Z_', frame,'.fits'), ext=1)
    imageRGB=list(R=Zimage$imDat, G=rimage$imDat, B=gimage$imDat)
    cat('images read in \n')
  
  
  
    SegimCat$segstats$starmask=SegimCat$segstats$xmax*0
    GroupCoordinates=Rwcs_p2s(SegimCat$segstats$xmax, SegimCat$segstats$ymax, header=Zimage$hdr)
    SegimCat$segstats$RAmax=GroupCoordinates[,'RA']
    SegimCat$segstats$Decmax=GroupCoordinates[,'Dec']
  
    ddfra=c(min(SegimCat$segstats$RAmax),max(SegimCat$segstats$RAmax))
    ddfdec=c(min(SegimCat$segstats$Decmax),max(SegimCat$segstats$Decmax))
  
    ra=median(SegimCat$segstats$RAmax)
    dec=median(SegimCat$segstats$Decmax)
    gaia0=gaia[ra > ddfra[1]-1.0 & ra < ddfra[2]+1.0 & dec > ddfdec[1]-1.0 & dec < ddfdec[2]+1.0,]
    gaia0=subset(gaia0,gaia0$phot_g_mean_mag<22)
    radius=4
    #
    maskedregion=coordmatch(coordcompare=SegimCat$segstats[,c("RAmax","Decmax")],coordref=gaia0[,list(ra,dec)],rad=radius,
      inunitref="deg",inunitcompare="deg",radunit="asec")
    maskedobjs=unique(maskedregion$ID[maskedregion$ID>0])
    SegimCat$segstats[maskedobjs,"starmask"]=1
  
    # now doing the same thing for the group stats
    SegimCat$groupstats$starmask=SegimCat$groupstats$xmax*0
    SegimCat$groupstats$matchedGAIAstar_mag=SegimCat$groupstats$xmax*0

    GroupCoordinates=Rwcs_p2s(SegimCat$groupstats$xmax, SegimCat$groupstats$ymax, header=Zimage$hdr)
    SegimCat$groupstats$RAmax=GroupCoordinates[,'RA']
    SegimCat$groupstats$Decmax=GroupCoordinates[,'Dec']
    #
    maskedregion=coordmatch(coordcompare=SegimCat$groupstats[,c("RAmax","Decmax")],coordref=gaia0[,list(ra,dec)],rad=radius,
      inunitref="deg",inunitcompare="deg",radunit="asec")
    maskedobjs=unique(maskedregion$ID[maskedregion$ID>0])
    Sel = maskedregion$ID[match(maskedobjs, maskedregion$ID)]
    SegimCat$groupstats[maskedobjs,"starmask"]=1
    temp=data.table(array1=maskedregion$ID[maskedregion$ID>0], array2=gaia0$phot_g_mean_mag[maskedregion$ID>0])
    MatchedGAIAstars=temp[,min(array2),by=array1]$V1
    SegimCat$groupstats[maskedobjs,"matchedGAIAstar_mag"]=MatchedGAIAstars
    # SegimCat$groupstats[maskedobjs,"matchedGAIAstar_mag"]=gaia0$phot_g_mean_mag[Sel]
    cat('GAIA stars identified \n')

    frameBorder=0 # change this value (was 850 for a full sq degree)
    ngroupThreshhold = 3
    magRange=c(9, 21.2)
  
    # I want to calculate the flux fracion of the leading segment of the total group
    SegimCat$groupstats$fluxFrac=SegimCat$groupstats$groupID*0
    SegimCat$groupstats$starAreaFrac=SegimCat$groupstats$groupID*0
    SegimCat$groupstats$nonStarMag=SegimCat$groupstats$groupID*0
    FluxScale=10^(-0.4*(0-8.9))
    SegimCat$groupstats$flux = SegimCat$groupstats$flux*FluxScale # rescaling the groupstats flux for the first batch of outputs, because this wasn't done internally to proFound (bug)
    for(ii in 1:length(SegimCat$groupstats$groupID[SegimCat$groupstats$starmask==1])){
      # want to calculate the fraction of starmask segments in total. 
      indices=1:length(SegimCat$groupstats$starmask)
      indices=indices[SegimCat$groupstats$starmask==1]
    
      StarFlux = sum(SegimCat$segstats$flux[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]] & SegimCat$segstats$starmask==1])
      StarArea = sum(SegimCat$segstats$N100[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]] & SegimCat$segstats$starmask==1])
      NonStarFlux=SegimCat$groupstats$flux[indices[ii]]-StarFlux
      NonStarMag=-2.5*log10(NonStarFlux)+8.9
  
      Fraction=StarFlux/(SegimCat$groupstats$flux[indices[ii]])
      AreaFraction=StarArea/(SegimCat$groupstats$N100[indices[ii]])
      
      SegimCat$groupstats$fluxFrac[indices[ii]]=Fraction
      SegimCat$groupstats$starAreaFrac[indices[ii]]=AreaFraction
      SegimCat$groupstats$nonStarMag[indices[ii]]=NonStarMag
  
    }
    cat('segment stellar stats determined \n')

    pro_r_detect = profoundProFound(rimage, segim=SegimCat$segim, static_photom=TRUE, verbose=TRUE, magzero=23.9, box=100)
    pro_Z_detect = profoundProFound(Zimage, segim=SegimCat$segim, static_photom=TRUE, verbose=TRUE, magzero=23.9, box=100)

    segstats = cbind(SegimCat$segstats,col=pro_r_detect$segstats$mag-pro_Z_detect$segstats$mag)

    groups = profoundAutoMerge(SegimCat$segim_orig, segstats = segstats, spur_lim=4e-3, col_lim = c(0,0.8), Ncut=1)
    segim_fix = profoundSegimKeep(SegimCat$segim, segID_merge=groups$segID)

    
    # selecting the groups which correspond to the automatically merged segments. 
    fixedSegments=unlist(groups$segID)
    fixedGroups = unique(SegimCat$segstats$groupID[match(fixedSegments, SegimCat$segstats$segID)]) # a list of groups that have an atomatic flag in them. 
    
    # now determining the final sample for viewing. 
    Sel = (SegimCat$groupstats$groupID%in%fixedGroups & 
		SegimCat$groupstats$xmax>=frameBorder & SegimCat$groupstats$xmax<=(dim(SegimCat$image)[1]-frameBorder) & 
		SegimCat$groupstats$ymax>=frameBorder & SegimCat$groupstats$ymax<=(dim(SegimCat$image)[1]-frameBorder) & 
      (SegimCat$groupstats$starmask==0 | (SegimCat$groupstats$starmask==1 & SegimCat$groupstats$starAreaFrac<0.7 & SegimCat$groupstats$nonStarMag<21.5 & SegimCat$groupstats$matchedGAIAstar_mag>11) ) )
    cat('final fixing sample determined \n')
    
    Coordinates = list(x=SegimCat$groupstats$xmax[Sel], y=SegimCat$groupstats$ymax[Sel])
  

    # if the N100 value is larger than some amount, then instead of having the coordinate at the brightest segment, 
    # instead have the coordinate to be the centre of the group. 
    indices=1:length(Sel)
    indices=indices[Sel]
    groupN100=SegimCat$groupstats$N100[Sel]
    groupStarmask=SegimCat$groupstats$starmask[Sel]
    groupGAIAmag=SegimCat$groupstats$matchedGAIAstar_mag[Sel]

    x_extra=c()
    y_extra=c()
    index_extra=c()

    if(length(indices)>0){
    	for(ii in 1:length(indices)){
    	  if(groupN100[ii]>100000){ # iterate to figure out what an appropriate value is here
    	    # in this case, update the coordinate values so that we are focusing on the centre of the group
    	    segIDs_group=SegimCat$segstats$segID[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]]
	
    	    # first checking whether two thumbnails are required to cover the whole group
    	    Minx=min(SegimCat$segstats$xmax[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]])
    	    Maxx=max(SegimCat$segstats$xmax[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]])
    	    Miny=min(SegimCat$segstats$ymax[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]])
    	    Maxy=max(SegimCat$segstats$ymax[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]])
	
    	    if( ((Maxx-Minx)>999) | ((Maxy-Miny)>999)){
    	      # now I need to determine the coordinates of each half of the object
    	      print(paste0('adding extra coordinate for group ', SegimCat$groupstats$groupID[indices[ii]]))
	
    	      if((Maxx-Minx)>(Maxy-Miny)){ # if the group is longer in x, then divide the groups horizontally
    	        Median=median(SegimCat$segstats$xmax[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]])
	
    	        segIDs_group_1=SegimCat$segstats$segID[(SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]) & (SegimCat$segstats$xmax>=Median)]
    	        segIDs_group_2=SegimCat$segstats$segID[(SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]) & (SegimCat$segstats$xmax<Median)]
	
    	      }else{ # if the group is longer in y, then divide the groups vertically
    	        Median=median(SegimCat$segstats$ymax[SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]])
	
    	        segIDs_group_1=SegimCat$segstats$segID[(SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]) & (SegimCat$segstats$ymax>=Median)]
    	        segIDs_group_2=SegimCat$segstats$segID[(SegimCat$segstats$groupID==SegimCat$groupstats$groupID[indices[ii]]) & (SegimCat$segstats$ymax<Median)]
    	      }
	
    	      # determining the central coordinates for the first half of the group
    	      x_coords_1=segIDs_group_1*0
    	      y_coords_1=segIDs_group_1*0
    	      subgroupN100_1=0 # the total pixel size of this half of the group 
    	      for(jj in 1:length(segIDs_group_1)){
    	        x_coords_1[jj]=SegimCat$segstats$xmax[SegimCat$segstats$segID==segIDs_group_1[jj]]*SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group_1[jj]]
    	        y_coords_1[jj]=SegimCat$segstats$ymax[SegimCat$segstats$segID==segIDs_group_1[jj]]*SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group_1[jj]]
    	        subgroupN100_1=subgroupN100_1+SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group_1[jj]]
    	      }
	
    	      # now replacing the main entry with these coordinates
    	      Coordinates$x[ii] = sum(x_coords_1)/subgroupN100_1
    	      Coordinates$y[ii] = sum(y_coords_1)/subgroupN100_1
	
    	      # determining the central coordinates for the second half of the group
    	      x_coords_2=segIDs_group_2*0
    	      y_coords_2=segIDs_group_2*0
    	      subgroupN100_2=0 # the total pixel size of this half of the group 
    	      for(jj in 1:length(segIDs_group_2)){
    	        x_coords_2[jj]=SegimCat$segstats$xmax[SegimCat$segstats$segID==segIDs_group_2[jj]]*SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group_2[jj]]
    	        y_coords_2[jj]=SegimCat$segstats$ymax[SegimCat$segstats$segID==segIDs_group_2[jj]]*SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group_2[jj]]
    	        subgroupN100_2=subgroupN100_2+SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group_2[jj]]
    	      }
	
    	      # saving these coordinates to be added back into the list after checking all objects
    	      x_extra = c(x_extra, sum(x_coords_2)/subgroupN100_2)
    	      y_extra = c(y_extra, sum(y_coords_2)/subgroupN100_2)
    	      index_extra = c(index_extra, ii)
	
    	    }else{ # in the case that the group will fit into a single thumbnail, simply calculate the weighted centre of the group. 
    	      print(paste0('changing coordinates for group ', SegimCat$groupstats$groupID[indices[ii]]))
    	      # now deriving the N100 of each segment
    	      x_coords=segIDs_group*0
    	      y_coords=segIDs_group*0
    	      for(jj in 1:length(segIDs_group)){
    	        x_coords[jj]=SegimCat$segstats$xmax[SegimCat$segstats$segID==segIDs_group[jj]]*SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group[jj]]
    	        y_coords[jj]=SegimCat$segstats$ymax[SegimCat$segstats$segID==segIDs_group[jj]]*SegimCat$segstats$N100[SegimCat$segstats$segID==segIDs_group[jj]]
    	      }
    	      Coordinates$x[ii] = sum(x_coords)/groupN100[ii]
    	      Coordinates$y[ii] = sum(y_coords)/groupN100[ii]
    	    }
  	
    	    
    	  }
    	}
    

    	# If extra coordinates have been saved, add them back into the main list. 
    	if(length(index_extra)>0){
    	  for(ii in 1:length(index_extra)){
    	    extra_index=index_extra[ii]
    	    groupN100=c(groupN100[1:extra_index], groupN100[extra_index], groupN100[(extra_index+1):length(groupN100)]) # replicate value
    	    # groupStarmask=c(groupStarmask[1:extra_index], groupStarmask[extra_index], groupStarmask[(extra_index+1):length(groupStarmask)])# replicate value
    	    # groupGAIAmag=c(groupGAIAmag[1:extra_index], groupGAIAmag[extra_index], groupGAIAmag[(extra_index+1):length(groupGAIAmag)])# replicate value
    	    Coordinates$x=c(Coordinates$x[1:extra_index], x_extra[ii], Coordinates$x[(extra_index+1):length(Coordinates$x)]) # add value
    	    Coordinates$y=c(Coordinates$y[1:extra_index], y_extra[ii], Coordinates$y[(extra_index+1):length(Coordinates$y)]) # add value
    	  }
    	}
    	cat('new thumbnail coordinates determined \n')
    	
  	
  	
    	box=ceiling(3*sqrt(groupN100))
    	box = box + (1 - (box%%2))
    	box[box<200]=200
    	box[box>999]=999 
    	
    	cat('Starting image saving process \n')
  	
    	dir.create(paste0(OutputDir, frame))
  	
    	print(paste('number of files', length(Coordinates$x)))
    	for(ii in 1:length(Coordinates$x)){
	
    	  print(paste(ii,'of',length(Coordinates$x)))
  	
    	  imageR=magcutout(imageRGB$R, loc=c(Coordinates$x[ii],Coordinates$y[ii]), box=box[ii])$image
    	  imageG=magcutout(imageRGB$G, loc=c(Coordinates$x[ii],Coordinates$y[ii]), box=box[ii])$image
    	  imageB=magcutout(imageRGB$B, loc=c(Coordinates$x[ii],Coordinates$y[ii]), box=box[ii])$image
    	  segim_cutout=magcutout(SegimCat$segim, loc=c(Coordinates$x[ii],Coordinates$y[ii]), box=box[ii])$image
    	  
  	
    	  Suffix=paste0('_', ii, '.fits')
	
    	  imageR[is.na(imageR)] = 0
    	  imageG[is.na(imageG)] = 0
    	  imageB[is.na(imageB)] = 0
  	
    	  Rfits_write_image(segim_cutout, filename=paste0(OutputDir, frame, '/', frame, Suffix))
    	  Rfits_write_image(imageR, filename=paste0(OutputDir, frame, '/', frame, Suffix,'[compress]'), create_file=FALSE)
    	  Rfits_write_image(imageG, filename=paste0(OutputDir, frame, '/', frame, Suffix,'[compress]'), create_file=FALSE)
    	  Rfits_write_image(imageB, filename=paste0(OutputDir, frame, '/', frame, Suffix,'[compress]'), create_file=FALSE)
    	}
	
    	saveRDS(list(segID_merge=groups$segID), paste0(OutputDir, 'waves_segID_merge_auto_', frame, '.rds'), version=2)
  	
    	# # now zipping up the directory
    	# files2zip <- dir(paste0(frame), full.names = TRUE)
    	# zip(zipfile = paste0(frame, '.zip'), files = files2zip)
    }else{ # in this case, there are no objects that need fixing! (somehow)
      print('no objects flagged for autmatic fixing')
      saveRDS(list(segID_merge=c()), paste0(OutputDir, 'waves_segID_merge_auto_', frame, '.rds'), version=2)

    }
  
    rm(SegimCat)
    rm(rimage)
    rm(gimage)
    rm(Zimage)
  }else{
    print(paste0(frame, ' output exists or input does not exist'))
  }

}



