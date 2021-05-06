from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.nddata.utils import Cutout2D

def bin_image(file_path = None, ext = None, use_array = False, data_array = None, bin_size=None, bin_number=None, bin_origin = 'lower left', fig_dimensions = (6,6), show_image = True, is_segm =  False, segm_object = None, copy = True, ignore_borders = True, border = 100):
    
    '''
    Bins fits images into bins of certain sizes of a specific number of bins then plots the resulting bins.
    
    Parameters
    -----------
    file_path - str
        File path of fits file
    ext - int
        Extension of fits file that contains the data
    bin_size - tuple (int)
        (xsize, ysize) of bins
    bin_number - tuple (int)
        (x_number, y_number) of bins
    bin_origin - str 
        Either 'lower left', 'upper right', or 'center' to specify where the bins align
        
    Outputs
    -------
    data - arr
        Numpy array of original data
    cutouts - list
        List of Cutout2D objects that are copies of the original data
    '''
    
    if use_array == True:
        data = data_array
        
    else:
        #Read file and data from file
        file = fits.open(file_path)
        data = file[ext].data
    
    #Raise exception if both bin size and bin number are defined
    if (not bin_size == None) and (not bin_number == None):
        raise Exception('Choose only bin size of bin number, not both.')
        
    if ignore_borders == False:
        border = 0

    img_xsize = len(data[0])-border*2
    img_ysize = len(data)-border*2
    
#    print(bin_size)
    
    #Based on wanted bin size, round DOWN to the number of bins that is closest to that bin size but still fills the image (with border) well
    #Want to do this so the border around the image is even
    if not bin_size == None:
  
        num_xbins = math.floor(img_xsize/bin_size[0])
        num_ybins = math.floor(img_ysize/bin_size[1])
        xsize = img_xsize/num_xbins
        ysize = img_ysize/num_ybins
        
#        print('num_xbins, num_ybins: {}, {}'.format(num_xbins, num_ybins))
#        print('xsize, ysize: {}, {}'.format(xsize, ysize))
        
    #Find size of bins based on bin number
    if not bin_number == None:
        num_xbins = bin_number[0]
        num_ybins = bin_number[1]
        xsize = img_xsize/bin_number[0]
        ysize = img_ysize/bin_number[1]
        bin_origin = 'lower left'
        
#    if (xsize > border) | (ysize > border):
#        raise Exception('Border size ({} px) is too small or bin size ({:.2f} px by {:.2f} px) is too big.'.format(border, xsize, ysize))
    
    #First cutout aligns with lower left corner
    if bin_origin == 'lower left':
        cutouts = [] #Make list of cutouts
        xstart = xsize/2 #Start cutout so the left side of the cutout aligns with x=0
        for xbin_index in range(0, num_xbins):
            ystart = ysize/2 #Start cutout so the bottom of the cutout aligns with y=0
            for ybin_index in range(0, num_ybins):
                cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=copy, mode='partial') #Make cutout
                ystart += ysize #Redefine position of cutout for next cutout
                cutouts.append(cutout)
            xstart += xsize
       
    #First cutout aligns with upper right corner
    if bin_origin == 'upper right':
        cutouts = []
        xstart = img_xsize-xsize/2 #Start cutout so the right side of the cutout aligns with x=xmax
        for xbin_index in range(0, num_xbins):
            ystart = img_ysize-ysize/2 #Start cutout so the top side of the cutout aligns with y=ymax
            for ybin_index in range(0, num_ybins):
                cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=copy, mode='partial')
                ystart -= ysize
                cutouts.append(cutout)
            xstart -= xsize
            
    #Cutouts are centered on image
    if bin_origin == 'center':
        cutouts = []
        
        if ignore_borders == False:
            xstart = (img_xsize % xsize)/2 - xsize/2
        xrange = num_xbins+1
        yrange = num_ybins+1
        if ignore_borders == True:
            xstart = border
#            xrange = num_xbins-2
#            yrange = num_ybins-2
            
#        print(xstart)
        
        for xbin_index in range(0, xrange):
        
            if ignore_borders == False:
                ystart = (img_ysize % ysize)/2 - ysize/2
            if ignore_borders == True:
                ystart = border
            
            for ybin_index in range(0, yrange):
#                print(xstart, ystart)
                
                if ignore_borders == False:
                    cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=True, mode='partial')
                if ignore_borders == True:
                    cutout = Cutout2D(data, (xstart, ystart), (ysize, xsize), copy=True, mode='partial')
                ystart += ysize
                cutouts.append(cutout)
            xstart += xsize
          
    if show_image == True:
        if is_segm == False:
            plt.figure(figsize=fig_dimensions)
            maskimg2 = np.ravel(data)
            maskimg2 = maskimg2[maskimg2 < 0.5] 
            minv2 = np.percentile(data,5)
            maxv2 = np.percentile(data,95)
            plt.imshow(data, vmin=minv2, vmax=maxv2, cmap='Greys_r', origin='lower')  
            for c in cutouts:
                c.plot_on_original(color='orange')
            plt.show()
        if is_segm == True:
            plt.figure(figsize=fig_dimensions)
            plt.imshow(segm_object, cmap=segm_object.cmap(random_state=12345), origin='lower')  
            for c in cutouts:
                c.plot_on_original(color='white') 
            plt.show()
    
    return (xsize, ysize), cutouts
