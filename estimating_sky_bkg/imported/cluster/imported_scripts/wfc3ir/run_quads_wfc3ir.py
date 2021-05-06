# Calculate sky of each WFC3IR quadrant for all images

import pandas as pd
from glob import glob
from astropy.io import fits

from percentileclip_quads_wfc3ir import calculate_sky

def calc_sky_quads(flter, file_list, px_area_arr):

    tempdf = pd.DataFrame([])

    for i, file in enumerate(file_list):

        print(i, file)

        fullsky, fullrms, quad1, quad2, quad3, quad4 = calculate_sky(file, 1, show = False, px_area_map = px_area_arr)

        quad1_sky = quad1[0]
        quad1_rms = quad1[1]
        quad1_good = quad1[2]
        quad1_bad = quad1[3]

        quad2_sky = quad2[0]
        quad2_rms = quad2[1]
        quad2_good = quad2[2]
        quad2_bad = quad2[3]

        quad3_sky = quad3[0]
        quad3_rms = quad3[1]
        quad3_good = quad3[2]
        quad3_bad = quad3[3]

        quad4_sky = quad4[0]
        quad4_rms = quad4[1]
        quad4_good = quad4[2]
        quad4_bad = quad4[3]

        loopdf = pd.DataFrame([[file.split('/')[6].split('_')[0], fullsky, fullrms, quad1_sky, quad1_rms, quad1_good, quad1_bad,
                              quad2_sky, quad2_rms, quad2_good, quad2_bad, 
                              quad3_sky, quad3_rms, quad3_good, quad3_bad, 
                              quad4_sky, quad4_rms, quad4_good, quad4_bad]], 
                             columns = ['root', 'fullimg_sky', 'fullimg_rms','quad1_sky', 'quad1_rms', 'quad1_good', 'quad1_bad',
                              'quad2_sky', 'quad2_rms', 'quad2_good', 'quad2_bad', 
                              'quad3_sky', 'quad3_rms', 'quad3_good', 'quad3_bad', 
                              'quad4_sky', 'quad4_rms', 'quad4_good', 'quad4_bad'])

        tempdf = tempdf.append(loopdf)

        tempdf.to_csv('/store/skysurf/wfc3ir/{f}/wfc3ir_{f}_percentileclip_quad_output.csv'.format(f = flter),index = False)

if __name__ == '__main__':

     flter = input('Flter? ')
     file_list = glob('/store/skysurf/wfc3ir/{f}/new_data/*_new_flt.fits.gz'.format(f = flter))
#     px_area_arr = fits.open('ir_wfc3_map.fits')[1].data

     calc_sky_quads(flter, file_list, px_area_arr)

