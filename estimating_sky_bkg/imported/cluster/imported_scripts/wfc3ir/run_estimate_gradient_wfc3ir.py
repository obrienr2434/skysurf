import pandas as pd
import os
import numpy as np
from glob import glob
from estimate_gradient_wfc3ir import calculate_gradient
import warnings
warnings.simplefilter(action='ignore')

camera = input('HST Camera (all lowercase): ')
flter = input('Filter (all uppercase): ')

cont = input('Is '+flter+' the correct filter? (yes or no) ')

if cont == 'yes':

        # Run code on all files in database for a particular filter
        # new_data corresponds to the data modified by Tim Carleton (up the ramp corrections)
        files = glob('/store/skysurf/'+camera+'/'+flter+'/new_data/*')

        df = pd.DataFrame([])

        j = 0
        for file in files:
                print(j, file)
                j +=1 # Keep track of the file # you are on

                root = os.path.basename(file.split('.')[0]) # Root name of the file (includes "new_flt")

                max_bkg, min_bkg, gradient = calculate_gradient(file, 1, sigm = 3, has_DQarr = True)

                tempdf = pd.DataFrame([[root,max_bkg,min_bkg,gradient]], 
                        columns = ['root','max_bkg','min_bkg','gradient'])

                df = df.append(tempdf)
                df.to_csv('~/gradients/'+camera+'_'+flter+'_gradient_estimate.csv', index = False)
