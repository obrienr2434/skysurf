# Create inspection plots for a certain filter

import pandas as pd
import os
import numpy as np
from glob import glob
from astropy.io import fits
from percentileclip_wfc3ir import calculate_sky

camera = input('HST Camera (all lowercase): ')
flter = input('Filter (all uppercase): ')

cont = input('Is '+flter+' the correct filter? (yes or no) ')

if cont == 'yes':

        # Make directory to store inspection plots
        plots_path = '/store/skysurf/'+camera+'/'+flter+'/inspection_plots/percentile-clip/'
        if os.path.exists(plots_path) == False:
                print('Making path '+plots_path)
                os.makedirs(plots_path)

        # Run code on all files in database for a particular filter
        # new_data corresponds to the data modified by Tim Carleton (up the ramp corrections)
        files = glob('/store/skysurf/'+camera+'/'+flter+'/new_data/*')

        print('N files = {}'.format(len(files)))

        j = 0
        for file in files:
                print(j, file)
                j +=1 # Keep track of the file # you are on

                root = os.path.basename(file.split('.')[0]) # Root name of the file (includes "new_flt")

                savepath = '/store/skysurf/'+camera+'/'+flter+'/inspection_plots/percentile-clip/'+root+'_'+flter+'_plot.png' # Where to save plots

                calc_bkg, calc_rms, N_good, N_bad, N_badpx, N_tot, expected_rms,exptime = calculate_sky(file, 1, 
                        savediag = True, savepath = savepath, show = False, has_DQ = True)

                tempdf = pd.DataFrame([[root,calc_bkg,calc_rms,N_good,N_bad,N_badpx,N_tot,expected_rms,exptime]], 
                        columns = ['root','calc_bkg','calc_rms','N_good','N_bad','N_badpx','N_tot','expected_rms','exptime'])

                df = df.append(tempdf)
                df.to_csv('/store/skysurf/'+camera+'/'+flter+'/'+camera+'_'+flter+'_percentileclip_output.csv', index = False)
