# plottools.py
# ALS 2016/04/28

"""
Tools for plotting images and meausrements (ellipse, contours, etc)

Contains two functions for plotting: plot_img, overplot_ellipse_fromtab

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import astropy.units as u
from astropy.table import Table

import polytools

def make_plot_img_w_contours(fn_plot, img, contours):
    """ 
    making a pdf plot at path fn_plot that have img in color and contours and color bar labeled by colorlabel 
    """
    fig, ax = plot_img(img, colorlabel=img.unit.to_string())
    overplot_contours(ax, contours)
    fig.savefig(fn_plot)


def make_iso_visual_panel(fn, img_compo, img_map, contours1, contours3, z, pixsize, legend_suffix, name, title_compo, title_map, label_cbar):
    fig, (ax0, ax1, ax_cb) = get_figure_grids(n_panel=2)

    # plot composit
    ax_imshow(fig, ax0, img_compo, origin='upper', tocolorbar=False, tosetlim=True)

    nx, ny = img_compo.shape[:2]
    overplot_ruler(ax0, z, pixsize=pixsize, rlength_arcsec=10., nx=nx, ny=ny)

    ax0.text(5, 12, name, color='white', fontsize=12)
    ax0.text(nx-35, 12, '$z={}$'.format('%.2f'%z), color='white', fontsize=10)
    ax0.set_title(title_compo)
    ax0.title.set_position([.5, 1.03])

    # plot line map
    im = ax_imshow(fig, ax1, img_map, vmin=-1, vmax=8, origin='lower', tocolorbar=False, tosetlim=True)
    overplot_contours(ax1, contours3, lw=1.)
    overplot_contours(ax1, contours1, lw=0.2)

    make_legend_isophotes(ax1, lw=2, suffix=legend_suffix)

    ax1.set_title(title_map)
    ax1.title.set_position([.5, 1.03])

    # plot color bar
    cbar = fig.colorbar(im, cax=ax_cb, label=label_cbar, format='%i')
    ax_cb.set_aspect(20)

    # set ticks off
    for ax in [ax0, ax1]: 
        ax.axis('off')

    # saving
    fig.savefig(fn, format='pdf')
    plt.close()


def get_figure_grids(n_panel=3):
    """ start figure and make subplots, there is one more additional (small) axis for colorbar """
    plt.close('all')
    fig=plt.figure(figsize=(2*n_panel+2.5, 3.))

    gs = gridspec.GridSpec(1, n_panel+1, width_ratios=[1 for i in range(n_panel)]+[0.05], wspace=0.1, left=0.05, right=0.90)

    axes = [plt.subplot(gs[i]) for i in range(n_panel+1)]

    return fig, axes


def plot_img(img, vmin=-1, vmax=10, origin='lower', tocolorbar=True, tosetlim=True, figsize=(6, 4.5), colorlabel='$I\/[10^{-15}\/\mathrm{erg\/\/s^{-1}\/cm^{-2}\/arcsec^{-2}}]$'):
    """
    Plot an image with colorbar. No saving is done. The fig and ax is 
    returned to make more modifications to the image. 

    PARAMS: 
        img
    OUTPUT: 
        fig, ax
    """
    plt.close('all')
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, aspect='equal')
    ax_imshow(fig, ax, img, vmin=vmin, vmax=vmax, origin=origin, tocolorbar=tocolorbar, tosetlim=tosetlim, colorlabel=colorlabel)
    ax.axis('off')    
    
    return fig, ax


def ax_imshow(fig, ax, img, vmin=None, vmax=None, origin='lower', tocolorbar=True, tosetlim=True, colorlabel = '$I\/[10^{-15}\/\mathrm{erg\/\/s^{-1}\/cm^{-2}\/arcsec^{-2}}]$'):
    """ for a given "ax" plot image"""
    im = ax.imshow(img, interpolation='nearest', origin=origin, aspect='equal', cmap='viridis', vmin=vmin, vmax=vmax)

    if tocolorbar:
        cb = fig.colorbar(im, label=colorlabel, format='%.1f')
        cb.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=15))
        cb.ax.tick_params(labelsize=12) 

    if tosetlim:
        ax.set_xlim(-0.5, img.shape[0]-0.5)
        if origin=='lower':
            ax.set_ylim(-0.5, img.shape[1]-0.5)
        else: 
            ax.set_ylim(img.shape[1]-0.5, -0.5)
    return im


def overplot_contour(ax, contour, color='white', ls='-',lw=3, alpha=1., label='__nolabel__'):
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
    ax.plot(contour[:, 1], contour[:, 0], linewidth=lw, color=color,lw=lw, alpha=alpha, label=label)



def overplot_contours(ax, contours, color='white', lw=3, alpha=1., label='__nolabel__'):
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
    for contour in contours:
        if polytools.SignedPolygonArea(contour)>0:
            overplot_contour(ax, contour, color=color, ls='-', lw=lw, alpha=alpha, label=label)
        else:
            overplot_contour(ax, contour, color=color, ls='--', lw=lw-1, alpha=alpha)



def overplot_ruler(ax, z, pixsize=0.396, rlength_arcsec=10., nx=64, ny=64):
    """
    Params
    ------
    ax:
    kpc_per_arcsec
    pixsize=0.396
        in arcsec
    rulerlength_arcsec=10.
        in arcsec
    """
    # import
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # set up
    xmid = 0.5*nx+5 # x ending point of bar
    y = ny-10. # y loc of bar

    # conversion
    kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(z)
    rlength_pix = rlength_arcsec/pixsize
    rlength_kpc = (rlength_arcsec*kpc_per_arcsec).value

    #===== plotting
    ax.plot([xmid-rlength_pix, xmid], [y, y], color='white', lw=2)
    ax.plot([xmid-rlength_pix, xmid-rlength_pix], [y+1., y-1.], color='white', lw=2)
    ax.plot([xmid, xmid], [y+1., y-1.], color='white', lw=2)
    ax.text(xmid+4., y+1, '10" ('+'%.0f'%rlength_kpc+' kpc)', color='white', fontsize=12)



def make_legend_isophotes(ax, lw=2, suffix=''):
    """ add legend of isophotes """
    ax.plot(-1, -1, color='white', lw=lw, label='Isophote'+suffix)
    # ax.plot(-1, -1, color='red', lw=lw, label='Galaxy mask')

    lg = ax.legend(loc='lower right', prop={'size':10})

    frame = lg.get_frame()
    frame.set_facecolor('#000099')
    frame.set_edgecolor('None')
    for tx in lg.get_texts():
        tx.set_color("white")

# def overplot_ellipse(ax, ellipse_params, color='black', label=None, lw=3):
#     """ overplot a ellipse on subplot ax. no units interpretation. 

#     Parameters
#     ----------
#     ax: matplotlib AxesSubplot object
#         to specify where to plot on

#     ellipse_params: list [xc, yc, a, b, theta]
#         xc: centroid x position
#         yc: centroid y position
#         a:  semi-major axis
#         b:  semi-minor axis  (all in pix)
#         theta: orientation of semi-major axis y of x in degrees

#     color: string 
#     """
#     from matplotlib.patches import Ellipse    

#     [xc, yc, a, b, theta]=ellipse_params

#     # define ellipse
#     kwrg = {'facecolor':'none', 'edgecolor':color, 'alpha':1, 'linewidth':lw}
#     ellip = Ellipse(xy=[xc,yc], width=2.*a, height=2.*b, angle=theta, **kwrg)
    
#     # overplotting and ellipse
#     ax.plot(xc, yc, marker='x', ms=5, mew=2, color=color)
#     ax.plot(0, 0, marker='', lw=lw, color=color, label=label)
#     ax.add_artist(ellip)




# def overplot_ellipse_fromtab(ax, img, tab_ellipse, color='black', pixelsize=0.396, useunits=True, label=None, lw=3):
#     """
#     PURPOSE: 
#         plot the corresponding ellipse on top of the image given ellipse 
#         measurements. a and b are treated as semi_major 
#         and semi_minor.  (diameter = 2.*a)

#     Parameters
#     ----------
#     ax: obj
#         the matplotlib subplot to overplot on

#     img: 2d np array 

#     tab_ellipse: astropy Table
#         table containing momoent measurements with columns including: 
#         'xc','yc','a','b','theta'
#         the quantities can have units. theta is in degrees. 

#     color='black': string

#     pixelsize=0.396: float
#         the pixelsize in arcsec if useunits is set to True

#     useunits=True: bool
#         If true, it assuems 'xc','yc','a','b' are in 
#         units of arcsec, and use pixelsize as specified. 
#     """

#     # unit conversion from arcsec to pix
#     if useunits:
#         for col in ['xc','yc','a','b']: # sanity check
#             if tab_ellipse[col].unit != u.Unit('arcsec'): 
#                 raise NameError('check unit')
#     else:
#         pixelsize=1.
#     xc = tab_ellipse['xc'][0]/pixelsize
#     yc = tab_ellipse['yc'][0]/pixelsize
#     a = tab_ellipse['a'][0]/pixelsize
#     b = tab_ellipse['b'][0]/pixelsize
#     theta = tab_ellipse['theta'][0]

#     ellipse_params=[xc, yc, a, b, theta]
#     overplot_ellipse(ax, ellipse_params, color=color, label=label, lw=lw)



# def overplot_axes(ax, params, color='black'):
#     """ overplot the two crossing perpendicular axes

#     Parameters
#     ----------
#     ax: matplotlib AxesSubplot object
#         to specify where to plot on

#     params: list [xc, yc, a, b, theta]
#         xc: centroid x position
#         yc: centroid y position
#         a:  major axis
#         b:  minor axis  (all in pix)
#         theta: orientation of major axis y of x in degrees

#     color: string 
#     """
#     from matplotlib.patches import Ellipse    

#     [xc, yc, a, b, theta]=params

#     va=0.5*a*np.array([np.cos(np.radians(theta)),np.sin(np.radians(theta))])
#     vb=0.5*b*np.array([-np.sin(np.radians(theta)),np.cos(np.radians(theta))])


#     # kwrg = {'color':color, 'lw':2}

#     ax.plot([xc-va[0],xc+va[0]],[yc-va[1],yc+va[1]],linewidth=2.0,color=color)#,*kwrg)
#     ax.plot([xc-vb[0],xc+vb[0]],[yc-vb[1],yc+vb[1]],linewidth=2.0,color=color)#,*kwrg)
#     ax.plot(xc,yc,marker='x',ms=5,mew=2,color=color)




# def overplot_feret_fromtab(ax, img, tab_feret, tab_xyc=None, color='black', pixelsize=0.396, useunits=True):
#     """
#     PURPOSE: 
#         plot the max feret distance and the max90 feret distance

#     Parameters
#     ----------
#     ax: obj
#         the matplotlib subplot to overplot on

#     img: 2d np array 

#     tab_feret: astropy Table
#         table containing columns including: 
#         'feretmax', 'theta_feretmax', 'feret90'
#         the quantities can have units. theta is in degrees. 

#     tab_xyc: astropy Table
#         table containing columns 'xc','yc'

#     color='black': string

#     pixelsize=0.396: float
#         the pixelsize in arcsec if useunits is set to True

#     useunits=True: bool
#         If true, it assuems distances are in 
#         units of arcsec, and use pixelsize as specified. 
#     """

#     # unit conversion from arcsec to pix
#     if useunits:
#         for col in ['feretmax','feret90']: # sanity check
#             if tab_feret[col].unit != u.Unit('arcsec'): 
#                 raise NameError('check unit')
#         for col in ['xc','yc',]: # sanity check
#             if tab_xyc is not None:
#                 if tab_xyc[col].unit != u.Unit('arcsec'): 
#                     raise NameError('check unit')
#     else:
#         pixelsize=1.

#     if tab_xyc is None:
#         xc, yc = standards.get_img_xycenter(img)        
#         # xc = 0.5*(img.shape[0]-1)
#         # yc = 0.5*(img.shape[0]-1)
#     else:
#         xc = tab_xyc['xc'][0]/pixelsize
#         yc = tab_xyc['yc'][0]/pixelsize
#     a=tab_feret['feretmax'][0]/pixelsize
#     b=tab_feret['feret90'][0]/pixelsize
#     theta = tab_feret['theta_feretmax'][0]

#     ellipse_params=[xc, yc, a, b, theta]
#     overplot_axes(ax, ellipse_params, color=color)


# def overplot_feretDR_frommdict(ax, dictparam, xc, yc, color='black'):
#     """
#     PURPOSE: 
#         plot the feretmaxD and feretmaxR
#     """

#     # unit conversion from arcsec to pix

#     vd=0.5*dictparam['dferetmax']*np.array([np.cos(np.radians(dictparam['theta_dferetmax'])),np.sin(np.radians(dictparam['theta_dferetmax']))])
#     vr=dictparam['rferetmax']*np.array([np.cos(np.radians(dictparam['theta_rferetmax'])),np.sin(np.radians(dictparam['theta_rferetmax']))])

#     ax.plot([xc-vd[0],xc+vd[0]],[yc-vd[1],yc+vd[1]],linewidth=4.0,color=color)
#     ax.plot([xc,xc+vr[0]],[yc,yc+vr[1]],linewidth=2.0,color=color)
#     ax.plot(xc,yc,marker='x',ms=5,mew=2,color=color)

