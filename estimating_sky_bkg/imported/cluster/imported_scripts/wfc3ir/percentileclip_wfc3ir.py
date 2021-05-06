import numpy as np
import os
from astropy.io import fits
from make_diagnostic import make_plots
from measureskyregion import measureskybin
from split_image import bin_image

def calculate_sky(filepath, ext, savediag = False, savepath = None, has_DQ = True,
RN_keyword = ['READNSEA','READNSEB','READNSEC','READNSED'], data_only = False):

	'''
	Params
	------
	filepath - str (fits)
		Path to fits file to be calculated
	ext - int
		Extension of SCI array
	savediag - bool
		Whether or not to save diagnostic/ inspection plot
	savepath - str
		Where to save diagnostic/ inspection plot
	has_DQ - bool
		Whether or not fits has a DQ extension
	RN_keyword - list of strings
		Header names for the readnoise for each amp on WFC3IR
	px_are_map - arr
		Pixel area map (multiplied by data)
	'''

	with fits.open(filepath) as hdu:

		data = hdu[ext].data
		data_wbadpix = np.copy(data)

		if has_DQ == True:
			data[hdu['DQ'].data != 0] = float('nan')

		shape = np.shape(data)
		if shape[0] == 1014:
			binpx = 39
			show_image = False
			# px_area_map_adj = px_area_map
		if (shape[0] == 512):
			binpx = 32
			# px_area_map_adj = px_area_map[251:763, 251:763]
		if (shape[0] == 256):
			binpx = 32
			# px_area_map_adj = px_area_map[379:635, 379:635]
		if (shape[0] == 128):
			binpx = 32
			# px_area_map_adj = px_area_map[443:571, 443:571]

		data = data*px_area_map_adj
		cutout_shape, cutouts = bin_image(use_array=True, data_array=data, bin_size=(binpx,binpx),
			bin_origin = 'lower left', show_image = False, ignore_borders = True, border = 0)
		
		#print(cutout_shape)

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

		#Good indices are where bkg array ISNT a NaN
		goodind = np.where(~np.isnan(bkg_arr))[0].tolist()

		#Define array without the nans
		bkg_arr_nonans = bkg_arr[goodind]
		calc_bkg = np.nanpercentile(bkg_arr_nonans,5) #Calc sky is 5th percentile of arr w/out nans
		calc_rms = np.nanmean(rms_arr) #Calc rms is just mean of all RMS values

		#Get INDICIES of lowest 5% of regions :)
		N_good = len(goodind) #Number of good regions
		N_5perc_of_good = int(N_good*0.05) #Number of lowest 5% of good regions
		lowest5perc = np.argsort(bkg_arr)[:N_5perc_of_good] #Indices of lowest 5% of good regions

		header = hdu[0].header
		targname = header['TARGNAME']

		if data_only == False:
			make_plots(data_wbadpix, cutouts, lowest5perc, badind, calc_bkg, calc_rms, badpx = badpxregs, calc_bkg = calc_bkg, calc_rms = calc_rms,
				gradient = None, title = os.path.basename(filepath), save = savediag, 
				savepath = savepath, show = False, targname = targname)

		exptime = header['EXPTIME']

		RN_list = []
		for kw in RN_keyword:
			rn = header[kw]
			RN_list.append(rn)

		RN = np.mean(RN_list)
		expected_rms = np.sqrt(RN**2+calc_bkg*exptime)/exptime

		hdu.close()

	return calc_bkg, calc_rms, len(lowest5perc), len(badind), len(badpxregs), len(cutouts), expected_rms, exptime