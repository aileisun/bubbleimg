# visual.py

"""
 make panel plots of composite and [OIII] line map, PSF residual
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.ndimage as simg
import pickle
import matplotlib as mpl

import astropy.table as at
from astropy.io import fits

from measureimg import plottools

# to embed fonts in pdf
import matplotlib
# matplotlib.rc('pdf', fonttype=42)

# dir_obj = '/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/mullaney/allT2_restricted/batch_rz/SDSSJ1352+6541/'

def batch_visual(dir_batch, dir_out, list_obj):
    """ make visual maps for obj in list_obj, data is from dir_bath/objname/, and output is under dir_out/objname_maps.pdf"""
    for objname in list_obj:
        dir_obj = dir_batch+objname+'/'
        filepath_out = dir_out+objname+'_maps.pdf'
        dir_makevisualmaps(dir_obj, filepath_out=filepath_out, objname=objname)


def batch_smallvisual(dir_batch, dir_out, list_obj):
    """ make visual maps for obj in list_obj, data is from dir_bath/objname/, and output is under dir_out/objname_maps.pdf"""
    for objname in list_obj:
        dir_obj = dir_batch+objname+'/'
        filepath_out = dir_out+objname+'_smallmaps.pdf'
        dir_makesmallvisualmaps(dir_obj, filepath_out=filepath_out, objname=objname)


def dir_makevisualmaps(dir_obj, filepath_out=None, objname=''):
    """ make panel representations of object """
    # set up 
    if filepath_out is None:
        filepath_out = dir_obj+'maps.pdf'

    # reading images
    img_gri = simg.imread(dir_obj+'HumVI_gri.png')
    img_lmap = fits.getdata(dir_obj+'stamp-lOIII5008_I_norm_denoised.fits')
    img_psfr = fits.getdata(dir_obj+'stamp-lOIII5008_I_norm_psfresidual_denoised.fits')

    # reading contours
    # ctr_lmap = read_pickle(dir_obj+'contours_blobctr.pkl')['contours']
    cdi_lmap = read_pickle(dir_obj+'contours_blobctr.pkl')
    isothreshold = cdi_lmap['isothreshold']
    ctr_lmap = cdi_lmap['contours']
    ctr_psfr = read_pickle(dir_obj+'contours_blob_psfresidctr.pkl')['contours']
    ctr_galm = read_pickle(dir_obj+'contours_galmaskctr.pkl')['contours']

    # read z
    z=at.Table.read(dir_obj+'xid.csv')['z'][0]

    # ==== plotting
    mpl.rcParams.update(mpl.rcParamsDefault)

    # setting 
    lw = 1.
    vmin = isothreshold*-0.5
    vmax = isothreshold*5.

    # makeing grid
    fig, ax0, ax1, ax2, ax3 = make_visualgrids()

    for ax in [ax0, ax1, ax2]: 
        ax.axis('off')    

    # plot gri png
    ax0.set_title('$\mathrm{Composite}\/-\/g\/r\/i}$')
    ax0.title.set_position([.5, 1.04])

    plottools.ax_imshow(fig, ax0, img_gri, origin='upper', tocolorbar=False, tosetlim=True)
    ax0.text(5, 8, objname, color='white')    
    overplot_ruler(ax0, z, pixsize=0.396, rlength_arcsec=10., nx=64, ny=64)

    # plot original line map
    ax1.set_title('$\mathrm{[OIII]\/\/Emission\/Line\/Map}$')
    ax1.title.set_position([.5, 1.04])
    plottools.ax_imshow(fig, ax1, img_lmap, vmin=vmin, vmax=vmax, tocolorbar=False, tosetlim=True)
    plottools.overplot_contours(ax1, ctr_galm, color='red', lw=lw)
    plottools.overplot_contours(ax1, ctr_lmap, color='white', lw=lw)

    # plot psf residual
    ax2.set_title('$\mathrm{[OIII]\/\/PSF\/Residual}$')
    ax2.title.set_position([.5, 1.04])
    im = plottools.ax_imshow(fig, ax2, img_psfr, vmin=vmin, vmax=vmax, tocolorbar=False, tosetlim=True)
    plottools.overplot_contours(ax2, ctr_psfr, color='white', lw=lw)
    make_legend_isophotes(ax2)

    # make color bar
    cbar = fig.colorbar(im, cax=ax3, label='$I\/[10^{-15}\/\mathrm{erg\/\/s^{-1}\/cm^{-2}\/arcsec^{-2}}]$', format='%i')
    ax3.set_aspect(20)
    # setting and saving
    # matplotlib.rc('pdf', fonttype=42)
    fig.savefig(filepath_out, format='pdf')


def make_visualgrids():
    """ start figure and make four subplots """
    plt.close('all')
    fig=plt.figure(figsize=(8.5, 3.))
    # partition panels
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.1, left=0.05, right=0.90)

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])

    return fig, ax0, ax1, ax2, ax3


def make_legend_isophotes(ax, lw=2):
    """ add legend of isophotes """
    ax.plot(-1, -1, color='white', lw=lw, label='Isophote')
    ax.plot(-1, -1, color='red', lw=lw, label='Galaxy mask')


    lg = ax.legend(loc='lower right', prop={'size':8})

    frame = lg.get_frame()
    frame.set_facecolor('#000099')
    frame.set_edgecolor('None')
    for tx in lg.get_texts():
        tx.set_color("white")


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
    xmid = 0.5*nx-2. # x ending point of bar
    y = ny-5. # y loc of bar

    # conversion
    kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(z)
    rlength_pix = rlength_arcsec/pixsize
    rlength_kpc = (rlength_arcsec*kpc_per_arcsec).value

    #===== plotting
    ax.plot([xmid-rlength_pix, xmid], [y, y], color='white', lw=2)
    ax.plot([xmid-rlength_pix, xmid-rlength_pix], [y+1., y-1.], color='white', lw=2)
    ax.plot([xmid, xmid], [y+1., y-1.], color='white', lw=2)
    ax.text(xmid+2., y+1, '10" ('+'%.0f'%rlength_kpc+' kpc)', color='white')


def read_pickle(filename):
    with open(filename, 'rb') as handle:
        result = pickle.load(handle)
    return result


def dir_makesmallvisualmaps(dir_obj, filepath_out=None, objname='', totitle=False):
    """ make panel representations of object """
    # set up 
    if filepath_out is None:
        filepath_out = dir_obj+'smallmaps.pdf'

    # reading images
    img_gri = simg.imread(dir_obj+'HumVI_gri.png')
    img_lmap = fits.getdata(dir_obj+'stamp-lOIII5008_I_norm.fits')

    # reading contours
    # ctr_lmap = read_pickle(dir_obj+'contours_blobctr.pkl')['contours']
    cdi_lmap = read_pickle(dir_obj+'contours_blobctr.pkl')
    isothreshold = cdi_lmap['isothreshold']

    # read z
    z=at.Table.read(dir_obj+'xid.csv')['z'][0]

    # ==== plotting
    mpl.rcParams.update(mpl.rcParamsDefault)

    # setting 
    lw = 1.
    vmin = isothreshold*-0.5
    vmax = isothreshold*5.

    # makeing grid
    fig, ax0, ax1 = make_smallvisualgrids()

    for ax in [ax0, ax1]: 
        ax.axis('off')    

    # plot gri png
    if totitle:
        ax0.set_title('$\mathrm{Composite}\/-\/g\/r\/i}$')
    plottools.ax_imshow(fig, ax0, img_gri, origin='upper', tocolorbar=False, tosetlim=True)
    ax0.text(5, 8, objname, color='white')    
    overplot_ruler(ax0, z, pixsize=0.396, rlength_arcsec=10., nx=64, ny=64)

    # plot original line map
    if totitle:
        ax1.set_title('$\mathrm{[OIII]\/\/Emission\/Line\/Map}$')
    plottools.ax_imshow(fig, ax1, img_lmap, vmin=vmin, vmax=vmax, tocolorbar=False, tosetlim=True)

    # setting and saving
    # matplotlib.rc('pdf', fonttype=42)
    
    fig.savefig(filepath_out, format='pdf')


def make_smallvisualgrids():
    """ start figure and make four subplots """
    plt.close('all')
    fig=plt.figure(figsize=(5.5, 3.))
    # partition panels
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1,], wspace=0.05, left=0.05, right=0.95)#, top=0.98, bottom=0.02)

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    return fig, ax0, ax1
