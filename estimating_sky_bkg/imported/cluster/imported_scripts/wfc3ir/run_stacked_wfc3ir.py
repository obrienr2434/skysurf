# Run percentileclip_stacked_wfc3ir.py

import pandas as pd
import os
import numpy as np
from astropy.io import fits
from percentileclip_stacked_wfc3ir import calculate_sky_stacked

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)

# Get sky and rms arrays from sky calculations to make stacked inspection plots
def get_stacked_data(flter, root_list, savepath = '~/python/', ramp_correct = True):

    '''
    Params
    ------
    flter - str
        WFC3IR filter
    root_list - list
        List of root names to return stacked data
    savepath - str
        Directory to save data
    ramp_correct - bool
        Whether or not to use up-the-ramp corrected images

    Outputs
    -------
    1) Saves each sky array, RMS array, and expected RMS array to data4stacking/ directory
    2) Saves a stacked SCI array of all the SCI data arrays
    '''
    
    df = pd.DataFrame([]) #Empty df to save sky and rms arrays
    
    data_stack = np.zeros((1014,1014)) #For stacking SCI data (for figure)
    
    #Loop through files in root_list
    #Enumerate starting at 1 so you can use variable "i" to stack data
    for i, root in enumerate(root_list, start = 1):
        
        print(i, root)
            
        if ramp_correct == True:
            filepath = '/store/skysurf/wfc3ir/{}/new_data/{}_new_flt.fits.gz'.format(flter, root)

        if ramp_correct == False:
            filepath = '/media/WD10TB_cl1/skysurf/wfc3ir/{}/data/{}_flt.fits.gz'.format(flter, root)

        skyarr, rmsarr, data, cutouts, exptime, readnoise = calculate_sky_stacked(filepath, 1)
        
        #Get expected rms of each region so that you can normalize RMS
        #So that final plot doesn't show sky bias
        expected_rms = np.sqrt(readnoise**2+skyarr*exptime)/exptime

        data_stack += data/i #stack SCI data (for figure)

        #make_dir(savepath)
        make_dir(os.path.join(savepath,'data4stacking'))
        make_dir(os.path.join(savepath,'data4stacking/{}_skyarr'.format(flter)))
        make_dir(os.path.join(savepath,'data4stacking/{}_rmsarr'.format(flter)))
        make_dir(os.path.join(savepath,'data4stacking/{}_expectedrmsarr'.format(flter)))

        np.savetxt(os.path.join(savepath,'data4stacking/{}_skyarr/{}_skyvals.txt'.format(flter, root)), skyarr)
        np.savetxt(os.path.join(savepath,'data4stacking/{}_rmsarr/{}_rmsvals.txt'.format(flter, root)), rmsarr)
        np.savetxt(os.path.join(savepath,'data4stacking/{}_expectedrmsarr/{}_expectedrmsarr.txt'.format(flter,root)), expected_rms)
    
    np.savetxt(os.path.join(savepath,'{}_datastack.txt'.format(flter)), data_stack) #Save data stack

if __name__ == '__main__':

    flter = input('Filter? ')

    root_list_path = '~/{}_quad_imageinfo.csv'.format(flter)
    savepath = 'quad_analysis/{}/'.format(flter)
    root_list = pd.read_csv(root_list_path)['root'].to_list()

    get_stacked_data(flter, root_list, savepath = savepath, ramp_correct = True)
