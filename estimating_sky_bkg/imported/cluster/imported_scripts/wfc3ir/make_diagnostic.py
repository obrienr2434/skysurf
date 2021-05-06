import matplotlib.pyplot as plt
import numpy as np

def make_plots(data, cutouts, goodind, badind, sky, rms, badpx = None, calc_bkg = None, calc_rms = None, gradient = None,
    title = None, save = False, savepath = None, show = False, targname = None):
    '''
    Params
    ------
    data - arr
        Full array containing SCI image
    cutoutskys - list
        List of 2DCutout objects
    goodind - list
        List of indices for "good" cutout regions
    badind - list
        List of indices for "bad" cutout regions
    '''
    
    reg_frac = len(badind)/len(cutouts)
    
    #Make sure path where image will save is specified if you are choosing to save the output
    if (save == True) & (savepath == None):
        raise Exception('Must specify savepath')

    fig, ax = plt.subplots(figsize = (7,7))
        
    if sky < 0:
        adjustdata = data+abs(np.percentile(data,1))
        minv = np.log10(sky+abs(np.percentile(data,1)))
        maxv =np.log10(sky+abs(np.percentile(data,1))+20*rms)
    elif rms < sky:
        minv = np.log10(sky-2*rms)
        maxv = np.log10(sky+20*rms)
        adjustdata = data+3*rms
    else:
        minv = np.log10(sky)
        maxv = np.log10(sky+20*rms)
        adjustdata = data+3*rms
        
    ax.imshow(np.log10(adjustdata), vmin=minv, vmax=maxv, cmap='Greys', origin='lower')
    
    #For each cutout object, plot it on the original image
    for ci, c in enumerate(cutouts):
        #Good regions are green
        if ci in goodind:
            c.plot_on_original(fill = False, color='green', label='Lowest 5% of good regions' if ci == goodind[0] else '', alpha = 0.5)
            c.plot_on_original(fill = True, facecolor='green', alpha = 0.09)
        #Bad regions are red
        if not np.shape(badind) == (0,): #Only run if there are bad regions
            if ci in badind:
                c.plot_on_original(fill = False, color='red', label='Bad regions' if ci == badind[0] else '', alpha = 0.5)
                c.plot_on_original(fill = True, facecolor='red', alpha = 0.09)
        if ci in badpx:
            c.plot_on_original(fill = False, color='purple', label='Regions with too many bad px' if ci == badpx[0] else '', alpha = 0.5)
            c.plot_on_original(fill = True, facecolor='purple', alpha = 0.09)
        #Neutral regions aren't colored
        if (not ci in goodind) and (not ci in badind) and (not ci in badpx):
            c.plot_on_original(color='white', alpha = 0.3) #White border
            
    #Add title to image
    if not title == None:
        plt.title(title)
        
#    print(reg_frac)

    #Add target name, calculated sky, calculated RMS, and calculated gradient
    if gradient == None:
        plt.text(0,-len(data)*0.125,'TARGNAME: {} \nsky: {:.3f} \nrms: {:.3f} \nF={:.3f}'.format(targname,calc_bkg,calc_rms,reg_frac))
    else:
        plt.text(0,-len(data)*0.125,'TARGNAME: {} \nsky: {:.2f} \nrms: {:.2f} \ngradient: {:.2f}'.format(targname,calc_bkg,calc_rms,gradient))

    #No tick marks
    ax.tick_params(labelbottom=False, labelleft = False)

    #Legend to indicate bad = red and good = green
    ax.legend(frameon = False, bbox_to_anchor=(1, -0.15), loc='lower right')
    
    if save == True:
        plt.savefig(savepath, bbox_inches = 'tight')
    if show == True:
        plt.show()
    if not show == True:
        plt.close()
