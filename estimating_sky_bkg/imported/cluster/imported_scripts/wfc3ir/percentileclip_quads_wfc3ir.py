import numpy as np
import os
from astropy.io import fits
from make_diagnostic import make_plots
from measureskyregion import measureskybin
from split_image import bin_image

# Split the Cutout2D objects so that they are grouped by quadrant (refer to wfc3ir instrument handbook)
def split_cutouts(arr):
    n_percol = int(np.sqrt(len(arr)))
    n_perrow = int(n_percol/2)
    
    n = n_percol
    columns = [arr[i:i + n] for i in range(0, len(arr), n)]
    
    split_data = []
    for column in columns:
        lst = column
        n = n_perrow
        rows = [lst[i:i + n] for i in range(0, len(lst), n)]
        split_data.append(rows)
    split_data = np.array(split_data)

    len_arr = len(arr)
#     print(len_arr)
    half_cols_ind = int(np.sqrt(len_arr)/2)
#     print(half_cols_ind)
    first_half_cols = split_data[0:half_cols_ind]
#     print(first_half_cols)
    sec_half_cols = split_data[half_cols_ind:int(2*half_cols_ind+1)]

    quad1 = first_half_cols[:,1] #top left
    quad2 = first_half_cols[:,0] #bottom left
    quad3 = sec_half_cols[:,0] #bottom right
    quad4 = sec_half_cols[:,1] #top right 

    return flatten(quad1), flatten(quad2), flatten(quad3), flatten(quad4)

def flatten(arr):
	return np.hstack(arr.flatten())

# Get sky, rms, and good/ bad region data
def sky_quad(bkg_arr, rms_arr, data_wbadpix, cutouts, bad_ind, targname, filepath, show):

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

#	make_plots(data_wbadpix, cutouts, lowest5perc, bad_ind, calc_bkg, calc_rms, badpx = bad_ind, calc_bkg = calc_bkg, calc_rms = calc_rms,
#		gradient = None, title = os.path.basename(filepath), save = False, 
#		savepath = None, show = show, targname = targname)

	return calc_bkg, calc_rms, len(lowest5perc), len(bad_ind)

def calculate_sky(filepath, ext, savediag = False, savepath = None, show = False, has_DQ = True,
RN_keyword = ['READNSEA','READNSEB','READNSEC','READNSED']):

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

		quad1_skys, quad2_skys, quad3_skys, quad4_skys = split_cutouts(bkg_arr)
		quad1_rms, quad2_rms, quad3_rms, quad4_rms = split_cutouts(rms_arr)
		quad1_cs, quad2_cs, quad3_cs, quad4_cs = split_cutouts(cutouts)

		header = hdu[0].header
		targname = header['TARGNAME']

		quad1_outputs = sky_quad(quad1_skys, quad1_rms, data_wbadpix, quad1_cs, np.where(np.isnan(quad1_skys))[0], targname, filepath, show)
		quad2_outputs = sky_quad(quad2_skys, quad2_rms, data_wbadpix, quad2_cs, np.where(np.isnan(quad2_skys))[0], targname, filepath, show)
		quad3_outputs = sky_quad(quad3_skys, quad3_rms, data_wbadpix, quad3_cs, np.where(np.isnan(quad3_skys))[0], targname, filepath, show)
		quad4_outputs = sky_quad(quad4_skys, quad4_rms, data_wbadpix, quad4_cs, np.where(np.isnan(quad4_skys))[0], targname, filepath, show)

		hdu.close()

		return calc_bkg, calc_rms, quad1_outputs, quad2_outputs, quad3_outputs, quad4_outputs
