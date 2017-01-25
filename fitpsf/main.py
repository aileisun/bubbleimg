# main.py

"""
for a dir_obj, find best fit psf to line map, subtract psf from line map.
"""
import os
import copy

import astropy.table as at
from astropy.io import fits

import loadpsf
import fitpsf


def dir_fit_psf(dir_obj, band='r', fileimg='stamp-lOIII5008_I_norm.fits', searchradius=5., fixb=True):
    """
    by default, fit r-band psf to map 'stamp-lOIII5008_I_norm.fits' The model has amplitute fixed to the maximum value within 5 pix of the img center, and the position is constrained to be within 5 pix as well. 
    """

    # define output filenames
    fileimg_base = os.path.splitext(fileimg)[0]

    filemodel = dir_obj+fileimg_base+'_psfmodel.fits'
    fileresid = dir_obj+fileimg_base+'_psfresidual.fits'
    fileparam = dir_obj+'psffit.csv'

    # read in
    psf = loadpsf.dir_load_psf(dir_obj, band=band)
    img = fits.getdata(dir_obj+fileimg)
    sigma = at.Table.read(dir_obj+'noiselevel.csv')['sigma'][0]

    # fit psf
    model, residual, params, chisq = fitpsf.fit_psf(img, psf, sigma,
                                        searchradius=searchradius, fixb=fixb)

    # prep output
    hdus_img = fits.open(dir_obj+fileimg)

    # output model
    hdus_model = copy.copy(hdus_img)
    hdus_model[0].data = model
    hdus_model[0].header['HISTORY'] = "Best fit PSF model"
    hdus_model.writeto(filemodel, clobber=True)

    # output residual
    hdus_resid = copy.copy(hdus_img)
    hdus_resid[0].data = residual
    hdus_resid[0].header['HISTORY'] = "Best fit PSF residual"
    hdus_resid.writeto(fileresid, clobber=True)

    # output params
    imgmax = img.max()
    modelmax = model.max()
    data = [[fileimg], [band], [imgmax], [modelmax], [fixb], [chisq],
            [params[0]], [params[1]], [params[2]]]
    names = ['fileimg', 'psfband', 'imgmax', 'modelmax', 'fixb', 'chisq', 'b', 'dx', 'dy']
    tabout = at.Table(data, names=names)
    tabout.write(fileparam)
