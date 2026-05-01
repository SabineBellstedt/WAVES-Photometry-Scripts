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
library(doParallel)
library(bit64)
CoreNumber = 10

inputargs=commandArgs(TRUE)
configFilename=as.character(inputargs[1])
config = yaml.load_file(configFilename)

InputDir = paste0(config$path$tile)
MeasureDir = paste0(config$path$measure, config$version$detect, config$version$measure, '/')
FibreMagDir = paste0(config$path$fibremag, config$version$detect, config$version$measure, '/')
PostprocessDir = paste0(config$path$postprocess, config$version$detect, config$version$measure, config$version$postprocess, '/')
CatsDir = paste0(config$path$catalogues, config$version$detect, config$version$measure, config$version$postprocess, config$version$catalogue, '/')
PlotDir = paste0(CatsDir, "plots/")

MeasureVersionSuffix = paste0('_', config$version$detect, config$version$measure)
PostprocessVersionSuffix = paste0('_', config$version$detect, config$version$measure, config$version$postprocess)

pixScale = 0.3
WAVESdepth = 21.1
fibreDiams = c(1.45)
zeropoint_image = list(g=0, r=0, Z=30)
zeropoint_original = list(g=0, r=0, Z=30) # as was used for the measure run

InputTargetCat=fread(paste0(config$path$reference, config$referencefiles$tilelist))

registerDoParallel(cores=CoreNumber)

foreach(j=(1:length(InputTargetCat$RA)))%dopar%{
# for(j in c(1)){

	ra=as.character((format(round(InputTargetCat$RA[j], 1), nsmall = 1)))
	dec=as.character((format(round(InputTargetCat$Dec[j], 1), nsmall = 1)))

	message(paste0('coordinates: ', ra, ' ', dec))
	OutputFilename = paste0(FibreMagDir, "waves_fibremags_", ra, '_', dec, MeasureVersionSuffix, '.parquet')
	if(!file.exists(OutputFilename)){
		MeasureFilename = paste0(MeasureDir,'waves_measured_', ra, '_', dec, MeasureVersionSuffix, '.rds')
		measure = readRDS(MeasureFilename)
	
		# need the sky-subtracted R-band image to do the fibre magnitude calculation
		image_g=Rfits_read_image(paste0(InputDir, 'VST/', config$version$stack, '/g_', ra, '_', dec, '_', config$version$stack, '.fits'))
		image_r=Rfits_read_image(paste0(InputDir, 'VST/', config$version$stack, '/r_', ra, '_', dec, '_', config$version$stack, '.fits'))
		image_Z=Rfits_read_image(paste0(InputDir, 'VISTA/', config$version$stack, '/Z_', ra, '_', dec, '_', config$version$stack, '.fits'))


		# convert pixels with exactly 0 flux to NA
		image_g$imDat[image_g$imDat==0L] = NA
		image_r$imDat[image_r$imDat==0L] = NA
		image_Z$imDat[image_Z$imDat==0L] = NA

		# computing the preferred cen and max coordinates to use. 

		tar_max = as.data.frame(list(segID = measure$cat_tot$segID, xcen=measure$cat_tot$xmax_rt, ycen=measure$cat_tot$ymax_rt))
		tar_cen = as.data.frame(list(segID = measure$cat_tot$segID, xcen=measure$cat_tot$xcen_rt, ycen=measure$cat_tot$ycen_rt))

		# now if sources are missing r-band information, then I want to replace the coordinate values with the delection band
		MissingSel = !is.finite(measure$cat_tot$flux_rt)
		tar_max$xcen[MissingSel] = measure$pro_detect$segstats$xmax[MissingSel]
		tar_max$ycen[MissingSel] = measure$pro_detect$segstats$ymax[MissingSel]

		tar_cen$xcen[MissingSel] = measure$pro_detect$segstats$xcen[MissingSel]
		tar_cen$ycen[MissingSel] = measure$pro_detect$segstats$ycen[MissingSel]

		coordinate_definition = rep('r', length(measure$cat_tot$segID))
		coordinate_definition[MissingSel] = 'detection'
	
		# now run the function 
		FibreMags_max_g = profoundAperPhot(image = image_g, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			keyvalues = NULL, 
			tar = tar_max,
			pixscale = pixScale, 
			magzero = zeropoint_image$g, 
			correction = TRUE, 
			# centype='max', 
			verbose = FALSE)
		FibreMags_cen_g = profoundAperPhot(image = image_g, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			keyvalues = NULL, 
			tar = tar_cen,
			pixscale = pixScale,
			magzero = zeropoint_image$g, 
			correction = TRUE, 
			# centype='mean', 
			verbose = FALSE)

		FibreMags_max_r = profoundAperPhot(image = image_r, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			keyvalues = NULL, 
			tar = tar_max,
			pixscale = pixScale, 
			magzero = zeropoint_image$r, 
			correction = TRUE, 
			# centype='max', 
			verbose = FALSE)
		FibreMags_cen_r = profoundAperPhot(image = image_r, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			keyvalues = NULL, 
			tar = tar_cen,
			pixscale = pixScale,
			magzero = zeropoint_image$r, 
			correction = TRUE, 
			# centype='mean', 
			verbose = FALSE)

		FibreMags_max_Z = profoundAperPhot(image = image_Z, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			keyvalues = NULL, 
			tar = tar_max,
			pixscale = pixScale, 
			magzero = zeropoint_image$Z, 
			correction = TRUE, 
			# centype='max', 
			verbose = FALSE)
		FibreMags_cen_Z = profoundAperPhot(image = image_Z, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			keyvalues = NULL, 
			tar = tar_cen,
			pixscale = pixScale,
			magzero = zeropoint_image$Z, 
			correction = TRUE, 
			# centype='mean', 
			verbose = FALSE)
	
		# subtract off the sky portion of the fibre mag based on measure$cat_tot$sky_mean_rt value in each segment
		# first compute the number of pixels in each fibre
		PixNumber_1 = (pi * (fibreDiams[1]/2)^2) / (pixScale^2)

		############ g-band ###############

		# now computing the mean sky flux in each source
		Fluxscale_sky = 10^(-0.4*(zeropoint_original$g-8.9)) # scale for converting sky vaues to Jansky
		Fluxscale_apPhot = 10^(-0.4*(zeropoint_image$g-8.9)) # scale for converting updated aperture photometry vaues to Jansky

		mean_sky_flux_1 = Fluxscale_sky * measure$cat_tot$sky_mean_gt * PixNumber_1
	
		# and now correcting all the columns
		FibreMags_max_g$flux_app_1_sub = Fluxscale_apPhot*FibreMags_max_g$flux_app_1 - mean_sky_flux_1 # in Jansky
		FibreMags_max_g$map_app_1_sub = 8.9 - 2.5*log10(FibreMags_max_g$flux_app_1_sub)

		FibreMags_cen_g$flux_app_1_sub = Fluxscale_apPhot*FibreMags_cen_g$flux_app_1 - mean_sky_flux_1  # in Jansky
		FibreMags_cen_g$map_app_1_sub = 8.9 - 2.5*log10(FibreMags_cen_g$flux_app_1_sub)

		############ r-band ###############

		# now computing the mean sky flux in each source
		Fluxscale_sky = 10^(-0.4*(zeropoint_original$r-8.9)) # scale for converting sky vaues to Jansky
		Fluxscale_apPhot = 10^(-0.4*(zeropoint_image$r-8.9)) # scale for converting updated aperture photometry vaues to Jansky

		mean_sky_flux_1 = Fluxscale_sky * measure$cat_tot$sky_mean_rt * PixNumber_1
	
		FibreMags_max_r$flux_app_1_sub = Fluxscale_apPhot*FibreMags_max_r$flux_app_1 - mean_sky_flux_1 # in Jansky
		FibreMags_max_r$map_app_1_sub = 8.9 - 2.5*log10(FibreMags_max_r$flux_app_1_sub)

		FibreMags_cen_r$flux_app_1_sub = Fluxscale_apPhot*FibreMags_cen_r$flux_app_1 - mean_sky_flux_1  # in Jansky
		FibreMags_cen_r$map_app_1_sub = 8.9 - 2.5*log10(FibreMags_cen_r$flux_app_1_sub)

		############ Z-band ###############

		# now computing the mean sky flux in each source
		Fluxscale_sky = 10^(-0.4*(zeropoint_original$Z-8.9)) # scale for converting sky vaues to Jansky
		Fluxscale_apPhot = 10^(-0.4*(zeropoint_image$Z-8.9)) # scale for converting updated aperture photometry vaues to Jansky

		mean_sky_flux_1 = Fluxscale_sky * measure$cat_tot$sky_mean_Zt * PixNumber_1

		FibreMags_max_Z$flux_app_1_sub = Fluxscale_apPhot*FibreMags_max_Z$flux_app_1 - mean_sky_flux_1 # in Jansky
		FibreMags_max_Z$map_app_1_sub = 8.9 - 2.5*log10(FibreMags_max_Z$flux_app_1_sub)

		FibreMags_cen_Z$flux_app_1_sub = Fluxscale_apPhot*FibreMags_cen_Z$flux_app_1 - mean_sky_flux_1  # in Jansky
		FibreMags_cen_Z$map_app_1_sub = 8.9 - 2.5*log10(FibreMags_cen_Z$flux_app_1_sub)
	
		# now including the postprocess output so that I can include the uberIDs
		PostprocessFilename = paste0(PostprocessDir,"waves_postprocessed_", ra, '_', dec, PostprocessVersionSuffix, '.parquet')
		postprocess = read_parquet(PostprocessFilename)

		# and now computing the mag errors in each tile
		FibreErrors_g = profoundAperRan(image = image_g, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			Nran = 1000,
			keyvalues = NULL, 
			pixscale = pixScale, 
			magzero = zeropoint_image$g, 
			correction = TRUE,
			fluxtype = 'Jansky', 
			verbose = FALSE)

		FibreErrors_r = profoundAperRan(image = image_r, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			Nran = 1000,
			keyvalues = NULL, 
			pixscale = pixScale, 
			magzero = zeropoint_image$r, 
			correction = TRUE,
			fluxtype = 'Jansky', 
			verbose = FALSE)

		FibreErrors_Z = profoundAperRan(image = image_Z, 
			segim = measure$pro_detect$dilated_segim, 
			app_diam = fibreDiams, 
			Nran = 1000,
			keyvalues = NULL, 
			pixscale = pixScale, 
			magzero = zeropoint_image$Z, 
			correction = TRUE,
			fluxtype = 'Jansky', 
			verbose = FALSE)
	
		Area = (pi * (fibreDiams[1]/2)^2)
	
		# now putting this all together into  single table to be saved
		output = data.table(uberID = as.integer64(postprocess$uberID), 
							segID = measure$cat_tot$segID, 
							xmax = tar_max$xcen, 
							ymax = tar_max$ycen, 
							xcen = tar_cen$xcen, 
							ycen = tar_cen$ycen, 
							coordinate_definition = coordinate_definition, 
							flux_cen_g_1p45 = FibreMags_cen_g$flux_app_1_sub,
							flux_cen_err_g_1p45 = rep(FibreErrors_g$errors[1], length(postprocess$uberID)),
							flux_cen_g_unsubtracted_1p45 = FibreMags_cen_g$flux_app_1,
							mag_cen_g_1p45 = FibreMags_cen_g$map_app_1_sub,
							SB_cen_g_1p45 = FibreMags_cen_g$map_app_1_sub+2.5*log10(Area),
							mag_cen_err_g_1p45 = (2.5/log(10))*(FibreErrors_g$errors[1]/FibreMags_cen_g$flux_app_1_sub),
							N_cen_g_1p45 = FibreMags_cen_g$N_app_1, 
							frac_cen_g_1p45 = FibreMags_cen_g$frac_app_1, 
							flux_max_g_1p45 = FibreMags_max_g$flux_app_1_sub,
							flux_max_err_g_1p45 = rep(FibreErrors_g$errors[1], length(postprocess$uberID)),
							flux_max_g_unsubtracted_1p45 = FibreMags_max_g$flux_app_1,
							mag_max_g_1p45 = FibreMags_max_g$map_app_1_sub,
							SB_max_g_1p45 = FibreMags_max_g$map_app_1_sub+2.5*log10(Area),
							mag_max_err_g_1p45 = (2.5/log(10))*(FibreErrors_g$errors[1]/FibreMags_max_g$flux_app_1_sub),
							N_max_g_1p45 = FibreMags_max_g$N_app_1, 
							frac_max_g_1p45 = FibreMags_max_g$frac_app_1, 
							flux_cen_r_1p45 = FibreMags_cen_r$flux_app_1_sub,
							flux_cen_err_r_1p45 = rep(FibreErrors_r$errors[1], length(postprocess$uberID)),
							flux_cen_r_unsubtracted_1p45 = FibreMags_cen_r$flux_app_1,
							mag_cen_r_1p45 = FibreMags_cen_r$map_app_1_sub,
							SB_cen_r_1p45 = FibreMags_cen_r$map_app_1_sub+2.5*log10(Area),
							mag_cen_err_r_1p45 = (2.5/log(10))*(FibreErrors_r$errors[1]/FibreMags_cen_r$flux_app_1_sub),
							N_cen_r_1p45 = FibreMags_cen_r$N_app_1, 
							frac_cen_r_1p45 = FibreMags_cen_r$frac_app_1, 
							flux_max_r_1p45 = FibreMags_max_r$flux_app_1_sub,
							flux_max_err_r_1p45 = rep(FibreErrors_r$errors[1], length(postprocess$uberID)),
							flux_max_r_unsubtracted_1p45 = FibreMags_max_r$flux_app_1,
							mag_max_r_1p45 = FibreMags_max_r$map_app_1_sub,
							SB_max_r_1p45 = FibreMags_max_r$map_app_1_sub+2.5*log10(Area),
							mag_max_err_r_1p45 = (2.5/log(10))*(FibreErrors_r$errors[1]/FibreMags_max_r$flux_app_1_sub),
							N_max_r_1p45 = FibreMags_max_r$N_app_1, 
							frac_max_r_1p45 = FibreMags_max_r$frac_app_1, 
							flux_cen_Z_1p45 = FibreMags_cen_Z$flux_app_1_sub,
							flux_cen_err_Z_1p45 = rep(FibreErrors_Z$errors[1], length(postprocess$uberID)),
							flux_cen_Z_unsubtracted_1p45 = FibreMags_cen_Z$flux_app_1,
							mag_cen_Z_1p45 = FibreMags_cen_Z$map_app_1_sub,
							SB_cen_Z_1p45 = FibreMags_cen_Z$map_app_1_sub+2.5*log10(Area),
							mag_cen_err_Z_1p45 = (2.5/log(10))*(FibreErrors_Z$errors[1]/FibreMags_cen_Z$flux_app_1_sub),
							N_cen_Z_1p45 = FibreMags_cen_Z$N_app_1, 
							frac_cen_Z_1p45 = FibreMags_cen_Z$frac_app_1, 
							flux_max_Z_1p45 = FibreMags_max_Z$flux_app_1_sub,
							flux_max_err_Z_1p45 = rep(FibreErrors_Z$errors[1], length(postprocess$uberID)),
							flux_max_Z_unsubtracted_1p45 = FibreMags_max_Z$flux_app_1,
							mag_max_Z_1p45 = FibreMags_max_Z$map_app_1_sub,
							SB_max_Z_1p45 = FibreMags_max_Z$map_app_1_sub+2.5*log10(Area),
							mag_max_err_Z_1p45 = (2.5/log(10))*(FibreErrors_Z$errors[1]/FibreMags_max_Z$flux_app_1_sub),
							N_max_Z_1p45 = FibreMags_max_Z$N_app_1,
							frac_max_Z_1p45 = FibreMags_max_Z$frac_app_1
							)
	
		setDF(output)
		# now saving the output as a .parquet file
		write_parquet(output, OutputFilename)
	}else{
		print('output already exists')
	}
}