# noiselevel.py
# ALS 2016/05/03

"""
estimate noise level on an image
"""
import os 
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import fits


def load_noiselevel(dir_obj,filename,update=True):
    """
    If dir_obj/'noiseleve.csv' exist and it has row with 
    filename=filename, then load the image sigma from the file. 
    Otherwise, calculate image noise level using getnoiselevel_gaussfit()

    Parameters
    ------
    dir_obj
    filename
    update=True
        if true, measure sigma and write/update 

    Returns
    ------
    sigma: (float)
    """
    # print "loading noiselevel"
    file_img=dir_obj+filename
    file_tab=dir_obj+'noiselevel.csv'
    if os.path.isfile(file_tab): 
        # noiselevel.csv exists
        tab=Table.read(file_tab)

        if (filename not in tab['filename']): 
            # measure sigma and write in table
            data=fits.open(file_img)[0].data
            sigma=getnoiselevel_gaussfit(data,dir_obj=dir_obj)
            if update:
                print "writing noise level"
                tabnew=Table([[filename],[sigma]],names=['filename','sigma'])
                tab=vstack([tab,tabnew])
                tab.write(file_tab,format='ascii.csv')
            return sigma

        elif (filename in tab['filename']) and update: 
            print "updating noise level"
            # measure sigma and update in table
            data=fits.open(file_img)[0].data
            sigma=getnoiselevel_gaussfit(data,dir_obj=dir_obj)
            # update
            tab['sigma'][tab['filename']==filename]=sigma
            tab.write(file_tab,format='ascii.csv')
            return sigma

        elif (filename in tab['filename']) and not update:
            # print "reading noise level from file"
            # read from file
            sigma=tab[tab['filename']==filename]['sigma'][0]
            if (type(sigma) is not float) and (type(sigma) is not np.float64): 
                raise NameError('Oh oh, sigma is not a float')
            return sigma

    else:
        # noiselevel.csv doesn't exist
        # measure sigma
        data=fits.open(file_img)[0].data
        sigma=getnoiselevel_gaussfit(data,dir_obj=dir_obj)
        print "writing noise level"    
        tab=Table([[filename],[sigma]],names=['filename','sigma'])
        tab.write(file_tab,format='ascii.csv')
        return sigma
            


def getnoiselevel_gaussfit(data,dir_obj,toplot=True):
    """
    return best fit gaussian noise level of the data, by fitting gaussian to
    the pixel value histogram. 

    Parameters:
    ------
    data: 2d np array
        image
    toplot: bool

    Returns
    ------
    sigma: float
    """
    from scipy.optimize import curve_fit

    # set up
    nb=1000 # number of bins

    histo,bin_edges=np.histogram(data.flatten(),bins=nb)

    bin_centers = bin_edges[0:-1]+np.diff(bin_edges)[0]/2.

    p0 = [np.max(histo), 0., np.std(data.flatten())] # A, x0, sigma
    coeff, var_matrix = curve_fit(gauss, bin_centers, histo, p0=p0)
    coeff[1:3]=coeff[1:3]
    sigma=np.abs(coeff[2])

    if toplot:
        # plotting
        plt.clf()
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(bin_centers,histo,label='data')
        ax.plot(bin_centers,gauss(bin_centers,*coeff),label='fit Gaussian sigma ='+'%.4e'%sigma)
        ax.legend()
        ax.set_xlabel('Intensity [img unit]')
        ax.set_ylabel('counts')
        fig.savefig(dir_obj+'noiselevel.pdf')
        plt.close()

    return sigma

def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


