# Imports
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import random
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.stats import mad_std
from photutils import make_source_mask
from photutils import MMMBackground
from photutils import SExtractorBackground
from photutils import MADStdBackgroundRMS
from scipy.stats import norm
import warnings
from split_image import bin_image
    
def get_rms(fits_data, sigm = 3, suppress_warnings = True):
    
    sigma_clip = SigmaClip(sigma=sigm)
    
    with warnings.catch_warnings():
        if suppress_warnings==True:
            warnings.simplefilter("ignore")

        bkgrms = MADStdBackgroundRMS(sigma_clip)
        bkgrms_value = bkgrms.calc_background_rms(fits_data)
        
        return bkgrms_value
        
def get_bkg_sext(fits_data, sigm = 3, suppress_warnings = True):

    sigma_clip = SigmaClip(sigma=sigm)

    with warnings.catch_warnings():
        if suppress_warnings==True:
            warnings.simplefilter("ignore")

        bkg = SExtractorBackground(sigma_clip)
        bkg_value = bkg.calc_background(fits_data)
            
        return bkg_value

    
def calculate_gradient(filepath, ext, sigm = 3, has_DQarr = False):
    '''
    Calculate sky-SB of fits image

    Parameters
    ----------
    filepath - str
        Path to fits file to calculate backgorund sky value
    ext - int
        Extension that contains data in fits file
    bins - tuple
        Number of bins in (x,y). If none, calculate_sky will not bin the image
    sigm - float
        Sigma value to be used for sigma-clipping
    method - str
        Photutils method to calculate sky background. Will throw error if 'MMM' or 'SExtractor' are not chosen because this script only implements these two methods
    has_DQarr - bool
        True if fits file contains a DQ array. Will assume the DQ array is the 3rd extension
    '''
    
    with fits.open(filepath) as hdu:
    
        flag = 0
    
        data = hdu[ext].data
        unmasked_data = np.copy(data)

        if has_DQarr == True:
            data[hdu['DQ'].data != 0] = float('nan')
        
        mask = make_source_mask(data, 1.5,
                    npixels=10, dilate_size=11)
        data[mask] = np.nan
        masked_arr = data
        
        shape = np.shape(data)
        if shape[0] == 1014:
            binpx = 39
        if (shape[0] == 512) | (shape[0] == 256):
            binpx = 32
        cutout_shape, cutouts = bin_image(use_array=True, data_array=data, bin_size=(binpx,binpx),
            bin_origin = 'lower left', show_image = False, ignore_borders = True, border = 0)
        
        # cutout_bkgTOT = []
        # cutout_rmsTOT = []
        og_cutout_bkgTOT = []
        og_cutout_rmsTOT = []
        
        bad_ind = []
        
        for ci, c in enumerate(cutouts):

            sky = get_bkg_sext(c.data, sigm = sigm)
            
            rms = get_rms(c.data, sigm = sigm)
            
            # cutout_bkgTOT.append(sky)
            # cutout_rmsTOT.append(rms)
            
            og_cutout_bkgTOT.append(sky)
            og_cutout_rmsTOT.append(rms)
    
        # for i, value in enumerate(cutout_bkgTOT):
        #     condition1 = ( value > min(np.array(cutout_bkgTOT)+np.array(cutout_rmsTOT))
        #     # condition2 = ( cutout_rmsTOT[i] > 2*np.nanmean(np.array(cutout_rmsTOT) ) )
        #     if condition1 :
        #         bad_ind.append(ci)
        
        # bkg_arr = np.array(cutout_bkgTOT)
        
        #Get max and min bkg for gradient calculation
        #Make array without nans (because nans go to "higher" end when sorting using np.sort)
        og_bkg_arr = np.array(og_cutout_bkgTOT)
        max_bkg = np.nanpercentile(og_bkg_arr,95)
        min_bkg = np.nanpercentile(og_bkg_arr,5)

        c1 = 0.738716326101489
        c0 = 2.4468018911443945
        gradient = 1/c1*(max_bkg-min_bkg)/min_bkg*100-c0

        hdu.close()

    return max_bkg, min_bkg, gradient
