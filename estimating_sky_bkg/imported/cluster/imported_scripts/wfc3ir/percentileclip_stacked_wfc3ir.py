# Calculates sky of each sub-region in an image and returns arrays of sky values and rms values

import numpy as np
import os
from astropy.io import fits
from measureskyregion import measureskybin
from split_image import bin_image

def calculate_sky_stacked(filepath, ext, has_DQ = True, RN_keyword = ['READNSEA','READNSEB','READNSEC','READNSED'], px_area_map = 1):

	'''
	Params
	------
	filepath - str (fits)
		Path to fits file to be calculated
	ext - int
		Extension of SCI array
	has_DQ - bool
		Whether or not fits has a DQ extension
	RN_keyword - list of strings
		Header names for the readnoise for each amp on WFC3IR
	px_are_map - arr
		Pixel area map (multiplied by data)

	Outputs
	-------
	bkg_arr - arr
		Sky values
	rms_arr - arr
		RMS values
	data - arr
		SCI data
	cutouts - list (Cutout2D objects)
	exptime - float
		Exposure time
	RN - float
		Readnoise
	'''

	with fits.open(filepath) as hdu:

		data = hdu[ext].data

		if has_DQ == True:
			data[hdu['DQ'].data != 0] = float('nan')

		# Only use images that are 1014x1014
		shape = np.shape(data)
		if shape[0] == 1014:
			binpx = 39
			cont = True
		if (shape[0] == 512) | (shape[0] == 256):
			cont = False
			print('!!! {} is not 1014 x 1014 px !!!'.format(filepath))
			hdu.close()

		if cont == True:

			cutout_shape, cutouts = bin_image(use_array=True, data_array=data, bin_size=(binpx,binpx),
				bin_origin = 'lower left', show_image = False, ignore_borders = True, border = 0)

			#Define list of sky and RMS values so it can easily be appended
			all_skys = []
			all_rms = []
			
			#Define list of bad pixel regions (will be based on DQ array)
			badpxregs = []

			for ci, c in enumerate(cutouts):
				sky, rms = measureskybin(c.data,axis=0)

				if np.count_nonzero(np.isnan(c.data)) > (cutout_shape[0]*cutout_shape[1])*0.2:
					sky = float('nan')
					rms = float('nan')
					badpxregs.append(ci)

				all_skys.append(sky)
				all_rms.append(rms)

			#Define all_skyprms as an array of each regions sky value plus its RMS value
			all_skyprms = np.array(all_skys)+np.array(all_rms)
			badind = []
			#Get list of bad indicies (that arent bad pixel values) and set these sky values to nan
			for sky_i, sky in enumerate(all_skys):
				if sky > np.nanmin(all_skyprms):
					badind.append(sky_i)
					all_skys[sky_i] = float('nan')
					all_rms[sky_i] = float('nan')

			#Define sky and RMS arrays after finished editing them
			bkg_arr = np.array(all_skys)
			rms_arr = np.array(all_rms)

			header = hdu[0].header
			exptime = header['EXPTIME']

			RN_list = []
			for kw in RN_keyword:
				rn = header[kw]
				RN_list.append(rn)

			RN = np.mean(rn)


			hdu.close()

			return bkg_arr, rms_arr, data, cutouts, exptime, RN
