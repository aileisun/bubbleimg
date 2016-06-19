# plottools.py
# ALS 2016/04/28

"""
Tools for plotting images and meausrements (ellipse, contours, etc)

Contains two functions for plotting: plot_img, overplot_ellipse_fromtab

"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table

from .. import standards


def plot_img(img,vmin=None,vmax=None):
    """
    Plot an image with colorbar. No saving is done. The fig and ax is 
    returned to make more modifications to the image. 

    PARAMS: 
        img
    OUTPUT: 
        fig, ax
    """

    plt.clf()
    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')
    cax=ax.imshow(img,interpolation='nearest',origin='lower left',aspect='equal',cmap='viridis',vmin=vmin,vmax=vmax)#,extent=[0,img.shape[0],0,img.shape[1]])
    
    fig.colorbar(cax)
    ax.set_xlim(-0.5,img.shape[0]-0.5)
    ax.set_ylim(-0.5,img.shape[1]-0.5)

    return fig, ax


def overplot_contour(ax, contour, color='black', ls='-',lw=3):
    """ 
    overplot a contour on subplot ax. no units interpretation. 

    Parameters
    ----------
    ax: matplotlib AxesSubplot object
        to specify where to plot on

    poly: (N*2) ndarray of double
        list of corner coordinates (r,c) of the contour 

    color: string 
    """
    ax.plot(contour[:, 1], contour[:, 0], linewidth=2, color=color,lw=lw)


def overplot_contours(ax, contours, color='black',lw=3):
    """ 
    overplot a set of contours on subplot ax. no units interpretation. 
    counter-clock wise contours will be solid lines while clock-wise contours 
    will be in dashed lines. 

    Parameters
    ----------
    ax: matplotlib AxesSubplot object
        to specify where to plot on

    poly: (N*2) ndarray of double
        list of corner coordinates (r,c) of the contour 

    color: string 
    """
    from polytools import SignedPolygonArea
    for contour in contours:
        if SignedPolygonArea(contour)>0:
            overplot_contour(ax, contour, color=color, ls='-',lw=lw)
        else:
            overplot_contour(ax, contour, color=color, ls='--',lw=lw-1)

def overplot_ellipse(ax, ellipse_params, color='black'):
    """ overplot a ellipse on subplot ax. no units interpretation. 

    Parameters
    ----------
    ax: matplotlib AxesSubplot object
        to specify where to plot on

    ellipse_params: list [xc, yc, a, b, theta]
        xc: centroid x position
        yc: centroid y position
        a:  semi-major axis
        b:  semi-minor axis  (all in pix)
        theta: orientation of semi-major axis y of x in degrees

    color: string 
    """
    from matplotlib.patches import Ellipse    

    [xc, yc, a, b, theta]=ellipse_params

    # define ellipse
    kwrg = {'facecolor':'none', 'edgecolor':color, 'alpha':1, 'linewidth':2}
    ellip = Ellipse(xy=[xc,yc], width=2.*a, height=2.*b, angle=theta, **kwrg)
    
    # overplotting and ellipse
    ax.plot(xc,yc,marker='x',ms=5,mew=2,color=color)
    ax.add_artist(ellip)




def overplot_ellipse_fromtab(ax,img,tab_ellipse,color='black',pixelsize=0.396,useunits=True):
    """
    PURPOSE: 
        plot the corresponding ellipse on top of the image given ellipse 
        measurements. a and b are treated as semi_major 
        and semi_minor.  (diameter = 2.*a)

    Parameters
    ----------
    ax: obj
        the matplotlib subplot to overplot on

    img: 2d np array 

    tab_ellipse: astropy Table
        table containing momoent measurements with columns including: 
        'xc','yc','a','b','theta'
        the quantities can have units. theta is in degrees. 

    color='black': string

    pixelsize=0.396: float
        the pixelsize in arcsec if useunits is set to True

    useunits=True: bool
        If true, it assuems 'xc','yc','a','b' are in 
        units of arcsec, and use pixelsize as specified. 
    """

    # unit conversion from arcsec to pix
    if useunits:
        for col in ['xc','yc','a','b']: # sanity check
            if tab_ellipse[col].unit != u.Unit('arcsec'): 
                raise NameError('check unit')
    else:
        pixelsize=1.
    xc = tab_ellipse['xc'][0]/pixelsize
    yc = tab_ellipse['yc'][0]/pixelsize
    a = tab_ellipse['a'][0]/pixelsize
    b = tab_ellipse['b'][0]/pixelsize
    theta = tab_ellipse['theta'][0]

    ellipse_params=[xc, yc, a, b, theta]
    overplot_ellipse(ax, ellipse_params, color=color)



def overplot_axes(ax, params, color='black'):
    """ overplot the two crossing perpendicular axes

    Parameters
    ----------
    ax: matplotlib AxesSubplot object
        to specify where to plot on

    params: list [xc, yc, a, b, theta]
        xc: centroid x position
        yc: centroid y position
        a:  major axis
        b:  minor axis  (all in pix)
        theta: orientation of major axis y of x in degrees

    color: string 
    """
    from matplotlib.patches import Ellipse    

    [xc, yc, a, b, theta]=params

    va=0.5*a*np.array([np.cos(np.radians(theta)),np.sin(np.radians(theta))])
    vb=0.5*b*np.array([-np.sin(np.radians(theta)),np.cos(np.radians(theta))])


    # kwrg = {'color':color, 'lw':2}

    ax.plot([xc-va[0],xc+va[0]],[yc-va[1],yc+va[1]],linewidth=2.0,color=color)#,*kwrg)
    ax.plot([xc-vb[0],xc+vb[0]],[yc-vb[1],yc+vb[1]],linewidth=2.0,color=color)#,*kwrg)
    ax.plot(xc,yc,marker='x',ms=5,mew=2,color=color)




def overplot_feret_fromtab(ax,img,tab_feret,tab_xyc=None,color='black',pixelsize=0.396,useunits=True):
    """
    PURPOSE: 
        plot the max feret distance and the max90 feret distance

    Parameters
    ----------
    ax: obj
        the matplotlib subplot to overplot on

    img: 2d np array 

    tab_feret: astropy Table
        table containing columns including: 
        'feretmax', 'theta_feretmax', 'feret90'
        the quantities can have units. theta is in degrees. 

    tab_xyc: astropy Table
        table containing columns 'xc','yc'

    color='black': string

    pixelsize=0.396: float
        the pixelsize in arcsec if useunits is set to True

    useunits=True: bool
        If true, it assuems distances are in 
        units of arcsec, and use pixelsize as specified. 
    """

    # unit conversion from arcsec to pix
    if useunits:
        for col in ['feretmax','feret90']: # sanity check
            if tab_feret[col].unit != u.Unit('arcsec'): 
                raise NameError('check unit')
        for col in ['xc','yc',]: # sanity check
            if tab_xyc is not None:
                if tab_xyc[col].unit != u.Unit('arcsec'): 
                    raise NameError('check unit')
    else:
        pixelsize=1.

    if tab_xyc is None:
        xc, yc = standards.get_img_xycenter(img)        
        # xc = 0.5*(img.shape[0]-1)
        # yc = 0.5*(img.shape[0]-1)
    else:
        xc = tab_xyc['xc'][0]/pixelsize
        yc = tab_xyc['yc'][0]/pixelsize
    a=tab_feret['feretmax'][0]/pixelsize
    b=tab_feret['feret90'][0]/pixelsize
    theta = tab_feret['theta_feretmax'][0]

    ellipse_params=[xc, yc, a, b, theta]
    overplot_axes(ax, ellipse_params, color=color)


def overplot_feretDR_frommdict(ax, dictparam, xc, yc, color='black'):
    """
    PURPOSE: 
        plot the feretmaxD and feretmaxR
    """

    # unit conversion from arcsec to pix

    vd=0.5*dictparam['dferetmax']*np.array([np.cos(np.radians(dictparam['theta_dferetmax'])),np.sin(np.radians(dictparam['theta_dferetmax']))])
    vr=dictparam['rferetmax']*np.array([np.cos(np.radians(dictparam['theta_rferetmax'])),np.sin(np.radians(dictparam['theta_rferetmax']))])

    ax.plot([xc-vd[0],xc+vd[0]],[yc-vd[1],yc+vd[1]],linewidth=4.0,color=color)
    ax.plot([xc,xc+vr[0]],[yc,yc+vr[1]],linewidth=2.0,color=color)
    ax.plot(xc,yc,marker='x',ms=5,mew=2,color=color)

