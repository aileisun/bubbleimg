# fitpsf.py

"""
fit psf to image 
"""
import scipy.ndimage.interpolation as interpolation
import numpy as np
import copy

import astropy.table as at
from astropy.io import fits

from .. import standards
import loadpsf

def fit_psf(img, psf, sigma, searchradius=None, fixb=False):
    """
    find best fit model to psf within search radius of centroid

    Params
    -------
    see fit_psf_params()

    Returns
    -----
    model: np 2d array
    residual: np 2d array
    params: [b, dx, dy]
        b: scaling of psf
        dx: shift of psf relative to img center
        dy: shift of psf relative to img center
    chisq: float
    """

    # fitting
    params, chisq = fit_psf_params(img, psf, sigma, 
                                   searchradius=searchradius, fixb=fixb)
    # calculate model and residual
    nx, ny = img.shape
    b, dx, dy = params
    model = shiftscalereshape_psf(psf, b, dx, dy, nx=nx, ny=ny)
    residual = img - model

    return model, residual, params, chisq


def fit_psf_params(img, psf, sigma, searchradius=None, fixb=False, verbose=False):
    """
    Find best fit psf model that fits the img within search radius. 
    Three params are considered: b - scale, dx - shift, dy - shift.

    The minimization starts from a init params of peak in the img.
    One can choose to fix psf amplitute to max of img by fixb=True. 

    Params
    -------
    img: 2d np array
    psf: 2d np array
    sigma: float
    searchradius: float in unit pix
        if set to float, consider psf only within a central radius of img. 
        default None - consider all image
    fixb=False: bool
        if true, b is fixed to max within search radius

    Returns 
    --------
    params_out: [b, dx, dy]
        b: scaling of psf
        dx: shift of psf relative to img center
        dy: shift of psf relative to img center
    chisq: float
    """
    from scipy.optimize import minimize

    # prep
    invsigmasq = 1./sigma**2
    nx, ny = img.shape
    xc, yc = standards.get_img_xycenter_fromnxny(nx, ny)

    # initializing params
    if searchradius is not None:
        # set search radius by zeroing img outside circle of radius r
        img_init = copy.copy(img)
        for i in range(nx):
            for j in range(ny):
                if (i-xc)**2+(j-yc)**2 > searchradius**2: 
                    img_init[i, j] = 0
        img_max=np.max(img_init)
        b0 = img_max/np.max(psf)
        x0, y0 = np.unravel_index(np.argmax(img_init), img_init.shape)-np.array([xc, yc])
    else: 
        img_max=np.max(img)
        b0 = np.max(img)/np.max(psf)
        x0, y0 = np.unravel_index(np.argmax(img), img.shape)-np.array([xc, yc])

    # fitting
    if img_max > 0:
        if fixb is False:
            results = minimize(chisq_bdxdy, x0=(b0, x0, y0), args=(img, psf,
                               invsigmasq, nx, ny), method='Powell')
            if verbose: print results
            params_out = results['x']
            chisq = results['fun']
        else: 
            results = minimize(chisq_dxdy, x0=(x0, y0), args=(b0, img, psf,
                               invsigmasq, nx, ny), method='Powell')
            if verbose: print results
            params_out = np.append(np.array([b0]), results['x'])
            chisq = results['fun']

    else: 
        print "skipping fit_psf_params as no signal in img"
        params_out = [0., 0., 0.]

    return params_out, chisq


#===== tool sets =====

def shiftscalereshape_psf(psf, b, dx, dy, nx=64, ny=64):
    """
    reshape psf (whatever dimension) to specified dimension (nx, ny), and
    shift the center to the center of the new image +dx +dy. The center of
    the  image is defined by standards.get_img_xycenter()
    Subpixel shift is done by third order spline interpolation.
    The extra paddding is filled with number 0.
    """

    # pad psf edges to match with img dimension
    dwx, dwy = np.array([nx, ny])-np.array(psf.shape)
    psfp = np.pad(psf, ((0, dwx), (0, dwy)),
                  mode='constant', constant_values=0.)
    # shift psf
    xc, yc = standards.get_img_xycenter_fromnxny(nx, ny)
    xc_psf, yc_psf = standards.get_img_xycenter(psf)
    dxc, dyc = xc-xc_psf, yc-yc_psf  # shift between ctrs
    psfc = interpolation.shift(psfp, shift=[dyc+dx, dxc+dy], order=3,
                               mode='constant', cval=0.0)
    return psfc*b


def chisq_bdxdy(bdxdy, img, psf, invsigmasq, nx=64, ny=64):
    """ return ((img - b*psf)/sigma)**2, where invsigmasq=1/sigma**2 """
    b, dx, dy = bdxdy
    psfc = shiftscalereshape_psf(psf, b, dx, dy, nx=nx, ny=ny)
    chisq = invsigmasq*np.sum((img - psfc)**2)
    return chisq


def chisq_dxdy(dxdy, b0, img, psf, invsigmasq, nx=64, ny=64):
    """ return ((img - b*psf)/sigma)**2, where invsigmasq=1/sigma**2 """
    dx, dy = dxdy
    b = b0
    bdxdy = [b, dx, dy]
    chisq = chisq_bdxdy(bdxdy, img, psf, invsigmasq, nx=nx, ny=ny)
    return chisq


#==== tests =====

def run_sim_for_test_fit_psf(psf, b=100., dx=5., dy=1., sigma=1., epsilons=(0.1, 0.5, 0.5), nx=64, ny=64):
    """
    run a simulation of creating shifted noised psf and retrieving
    the shifts, scales by fitting.

    Params
    ---------
    b: model = shifted psf * b
    epsilongs: error tolerance (frac error in b, abs error in dx, dy)
    """
    # make model
    params_in = np.array([b, dx, dy])
    img = noisedpsf(psf=psf, b=b, dx=dx, dy=dy, sigma=sigma,
                          nx=nx, ny=ny)

    params_out, __ = fit_psf_params(img, psf, sigma, searchradius=None)

    fracerrb = np.absolute(params_in[0]-params_out[0])/np.absolute(params_in[0])
    absoerr  = np.absolute(params_in-params_out)
    # pisgood = fracerr < epsilon
    pisgood = [fracerrb < epsilons[0], absoerr[1] < epsilons[1], 
               absoerr[2] < epsilons[2]]
    if np.all(pisgood):
        if verbose: print "Success"
        return True
    else:
        if verbose: print "Failed"
        return False


def main_run_sims_for_test_fit_psf(psf, brange = [10., 30., 50., 80., 100., 300., 1000.], dx=0., dy=0.):

    # set up
    nrun = 50
    kwargs = dict(sigma=1., epsilons=(0.1, 1., 1.), nx=64, ny=64)

    snrpeaks=np.zeros(len(brange))
    failedrates=np.zeros(len(brange))

    for i, b in enumerate(brange):
        snrpeak = psf.max()*b/sigma
        records = np.zeros(nrun, dtype=bool)
        for j in range(nrun):
            records[j] = run_sim(psf, b=b, dx=dx, dy=dy, **kwargs)
        failedrate = (nrun - records.sum())/float(nrun)

        snrpeaks[i] = snrpeak
        failedrates[i] = failedrate

        print "Params: b ", str(b), ", dx ", str(dx), ", dy ", str(dy)
        print "SNR Peak ", '%.1f'%snrpeak
        print "Ran ", str(nrun), " times"
        print "Failed rate ", failedrate*100., '%'

    tabout=at.Table([brange, snrpeaks, failedrates], 
                    names=['b', 'snrpeak', 'failedrate'])

    tabout.write('results_sim.csv', format='ascii.csv')


def noisedpsf(psf, b=100., dx=1., dy=1., sigma=1., nx=64, ny=64):
    """
    Return a model, a np array of shape (nx, ny), of a psf located at
    cetner +dx +dy, scaled by factor b and with noise of level sigma added.
    """
    psfc = shiftscalereshape_psf(psf, b, dx, dy, nx=64, ny=64)
    model = add_gaussnoise(img=psfc, noise_sigma=sigma)
    # clf()
    # imshow(model, interpolation='nearest', cmap='viridis')
    return model


def add_gaussnoise(img, noise_sigma):
    """ add iid gaussian noise on image"""
    from skimage.util import random_noise
    if noise_sigma == 0:
        return img
    else:
        return img+random_noise(np.zeros(img.shape), mode='gaussian', 
                                var=noise_sigma**2, clip=False)
