# main.py
# ALS 2016/04/27

"""
Interface with outer funcitons to handle all image measurements (moment, gaussfit) with specifics like units and pixel sizes. 

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle

import astropy.units as u
import astropy.table as at
from astropy.io import fits

import moment
import gaussfit
import iso
import plottools
import polytools
import shapelytools

# from bubblepy import standards
from .. import standards


def dir_doIsos(dir_obj, isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'),  isoareallimit=10, contrastr=0.1, update=True, toplot=True, toclean=False): 
    """
    do the following operations: 
        find/store isocut contour on denoised line map as linblob
        measure isoshape on the lineblob
        find/store isocut/contrastr contour on continuum map as galmask
        subtract galmask from lineblob as lineblobmasked
        measure isoshape on the lineblobmasked

    Notation:
        fc_: filename of contour dictionary
        c_: contour dictionary 
        blob: OIII line 
        gal: galaxy
        ctr: center (not including disconnected patches)
        mask: mask obtained by "contrastr" times higher isocut
        masked: applied mask subtraction
    """
    # ==== list of input filenames
    fileimg_line = 'stamp-lOIII5008_I_norm_denoised.fits'
    fileimg_line_psfresid = 'stamp-lOIII5008_I_norm_psfresidual_denoised.fits'
    fileimg_gal = 'stamp-conti-onOIIIscale_I_norm.fits'
    fileimg_psfmodel = 'stamp-lOIII5008_I_norm_psfmodel.fits'

    # ==== define contour names
    prefix_cdict='contours_'

    # subject contours
    fc_blob = prefix_cdict+'blob.pkl'
    fc_blobctr = prefix_cdict+'blobctr.pkl'
    fc_blob_psfresid = prefix_cdict+'blob_psfresid.pkl'
    fc_blob_psfresidctr = prefix_cdict+'blob_psfresidctr.pkl'
    fc_galctr = prefix_cdict+'galctr.pkl'
    fc_psf = prefix_cdict+'psf.pkl'

    # mask contours
    fc_galmask = prefix_cdict+'galmask.pkl'
    fc_galmaskctr = prefix_cdict+'galmaskctr.pkl'
    fc_psfmask = prefix_cdict+'psfmask.pkl'
    
    # masked contours
    fc_blob_psfresid_galm = prefix_cdict+'blob_psfresid_galmasked.pkl'
    fc_blob_psfresid_galm_psfm = prefix_cdict+'blob_psfresid_galmasked_psfmasked.pkl'
    fc_blob_psfresidctr_galm = prefix_cdict+'blob_psfresidctr_galmasked.pkl'
    fc_blob_psfresidctr_galm_psfm = prefix_cdict+'blob_psfresidctr_galmasked_psfmasked.pkl'

    # ==== clean files
    if toclean:
        fms=["measureimg_iso.csv", "measureimg_iso.ecsv",]
        fcs=[fc_blob, fc_blobctr, fc_blob_psfresid, fc_blob_psfresidctr, fc_galctr, fc_psf, fc_galmask, fc_galmaskctr, fc_psfmask, fc_blob_psfresid_galm, fc_blob_psfresid_galm_psfm, ]
        fplots=['measureimg_iso_'+fc.split('.')[0]+'.pdf' for fc in fcs]

        for filename in fms+fcs+fplots:
            if os.path.isfile(dir_obj+filename):
                print "deleting file "+filename
                os.remove(dir_obj+filename)
            else:
                print "skip deleting file "+filename

    # ==== make contours
    print "making contours"
    # subject contours
    c_blob = dir_makeContoursDict(dir_obj, fileimg=fileimg_line, filecontoursdict=fc_blob, isocut_rest=isocut_rest, isoareallimit=isoareallimit, update=update)

    c_blobctr = dir_makeContoursDict(dir_obj, fileimg=fileimg_line, filecontoursdict=fc_blobctr, isocut_rest=isocut_rest, isoareallimit=np.inf, update=update)

    c_blob_psfresid = dir_makeContoursDict(dir_obj, fileimg=fileimg_line_psfresid, filecontoursdict=fc_blob_psfresid, isocut_rest=isocut_rest, isoareallimit=isoareallimit, update=update)

    c_blob_psfresidctr = dir_makeContoursDict(dir_obj, fileimg=fileimg_line_psfresid, filecontoursdict=fc_blob_psfresidctr, isocut_rest=isocut_rest, isoareallimit=np.inf, update=update)

    c_galctr = dir_makeContoursDict(dir_obj, fileimg=fileimg_gal, filecontoursdict=fc_galctr, isocut_rest=isocut_rest, isoareallimit=np.inf, update=update)

    c_psf = dir_makeContoursDict(dir_obj, fileimg=fileimg_psfmodel, filecontoursdict=fc_psf, isocut_rest=isocut_rest, isoareallimit=0., update=update)

    # mask contours

    c_galmask = dir_makeContoursDict(dir_obj, fileimg=fileimg_gal, filecontoursdict=fc_galmask, isocut_rest=isocut_rest/contrastr, isoareallimit=0., update=update)

    c_galmaskctr = dir_makeContoursDict(dir_obj, fileimg=fileimg_gal, filecontoursdict=fc_galmaskctr, isocut_rest=isocut_rest/contrastr, isoareallimit=np.inf, update=update)

    c_psfmask = dir_makeContoursDict(dir_obj, fileimg=fileimg_psfmodel, filecontoursdict=fc_psfmask, isocut_rest=isocut_rest/contrastr, isoareallimit=0., update=update)

    # masked contours
    print "masking contours"

    c_blob_psfresid_galm = dir_makeDiffContoursDict(dir_obj, fc_blob_psfresid, fc_galmask, fc_blob_psfresid_galm, update=update)

    c_blob_psfresid_galm_psfm = dir_makeDiffContoursDict(dir_obj, fc_blob_psfresid_galm, fc_psfmask, fc_blob_psfresid_galm_psfm, update=update)

    c_blob_psfresidctr_galm = dir_makeDiffContoursDict(dir_obj, fc_blob_psfresidctr, fc_galmask, fc_blob_psfresidctr_galm, update=update)

    c_blob_psfresidctr_galm_psfm = dir_makeDiffContoursDict(dir_obj, fc_blob_psfresidctr_galm, fc_psfmask, fc_blob_psfresidctr_galm_psfm, update=update)

    # ==== make iso measurements
    print "make iso measurements"

    cdicts_tomsr=[c_blob, c_blobctr, c_blob_psfresid, c_blob_psfresidctr, c_galctr, c_psf, c_galmask, c_galmaskctr, c_psfmask, c_blob_psfresid_galm, c_blob_psfresid_galm_psfm, c_blob_psfresidctr_galm, c_blob_psfresidctr_galm_psfm, ]

    for cdict in cdicts_tomsr:
        dir_MeasureImgIso_fromContoursDict(dir_obj, cdict, towritetab=True, toplot=toplot, update=update)


def dir_makeDiffContoursDict(dir_obj, fcdict1, fcdict2, fcdict_diff, update=True):
    """
    write difference of contour dict to file
    """
    c1 = read_pickle(dir_obj+fcdict1)
    c2 = read_pickle(dir_obj+fcdict2)

    if not os.path.isfile(dir_obj+fcdict_diff) or update:
        cdiff = shapelytools.diffcontoursdict(c1, c2, fcdict_diff)
        write_pickle(cdiff, dir_obj+fcdict_diff)
    else: 
        cdiff = read_pickle(dir_obj+fcdict_diff)

    return cdiff


def read_pickle(filename):
    with open(filename, 'rb') as handle:
        result = pickle.load(handle)
    return result


def write_pickle(result, filename):
    with open(filename, 'wb') as handle:
        pickle.dump(result, handle)


def dir_makeContoursDict(dir_obj, fileimg='stamp-lOIII5008_I_norm.fits', filecontoursdict=None, isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'), isoareallimit=0, update=False):
    """
    return and write (optional) contours from img fits files

    Return 
    ------
    contours: list of poly
        the contours
    isothreshold: float
        the isothreshold in image units
    """
    # set up input filenames
    filein=dir_obj+fileimg
    # set up contour filenames
    tagisocut='_iso'+'%.0e'%isocut_rest.value
    if np.isinf(isoareallimit):
        tagalim='_alim'+'inf'
    else:
        tagalim='_alim'+str(int(isoareallimit))
    fileimg_base=os.path.splitext(fileimg)[0] 
    if filecontoursdict is None:
        filecontoursdict='contoursdict_'+fileimg_base+tagisocut+tagalim+'.pkl'
    fileout=dir_obj+filecontoursdict

    # operation
    if not os.path.isfile(fileout) or update:
        # read in images and isothresholds
        img=fits.getdata(filein)
        xc, yc = standards.get_img_xycenter(img)
        # xc, yc = 0.5*(np.array(img.shape)-1.) # not used here but stored
        isothreshold=dir_getIsothreshold(dir_obj, fileimg=fileimg, isocut_rest=isocut_rest,)
        # measurements
        contours=polytools.find_realisocontours(img, threshold=isothreshold, areallimit=isoareallimit)
        # constructing output dictionary
        contoursdict=dict(contours=contours, fileimg=fileimg, filecontoursdict=filecontoursdict, isocut_rest=isocut_rest, isoareallimit=isoareallimit, isothreshold=isothreshold, xc=xc, yc=yc, masked=False, mask='None')
        # writing
        with open(fileout, 'wb') as handle:
          pickle.dump(contoursdict, handle)

    else: 
        with open(fileout, 'rb') as handle:
          contoursdict=pickle.load(handle)

    return contoursdict


def dir_getIsothreshold(dir_obj, fileimg='stamp-lOIII5008_I_norm.fits', isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)')):
    """
    calculate isothreshold in native image units from isocut_rest

    Return 
    ------
    isothreshold: float
        isophotal threshold in native image unit
    """
    # read in imgunit
    filein=dir_obj+fileimg
    header=fits.getheader(filein)
    imgunit=u.Unit(header['BUNIT'])
    # read in z
    xid=at.Table.read(dir_obj+'xid.csv',format='ascii.csv',comment='#')
    z=xid['z'][0]
    # calcualte isothreshold
    isocut_obs=isocut_rest*(1.+z)**-4
    isothreshold=(isocut_obs/imgunit).to(u.dimensionless_unscaled).value
    return isothreshold


def dir_MeasureImgIso_fromContoursDict(dir_obj, contoursdict, towritetab=True, toplot=True, update=True):
    """
    Make iso measurements on the specified image in the directory, and write 
    the results in file dir_obj+'measureimg.csv'. 
    
    Parameters
    ------
    dir_obj: string
        path to directory

    filecontour
        ='contoursdict_stamp-lOIII5008_I_norm_denoised_iso3e-15_alim10.pkl'
        the contour dictionary file
    towritetab=True, toplot=True, update=False
    
    Returns
    ------
    None

    Write Output
    ------
    dir_obj+'measureimg.csv'
    """
    # set up filenames
    fileout=dir_obj+'measureimg_iso'
    fileplot=fileout+'_'+os.path.splitext(contoursdict['filecontoursdict'])[0] 

    # start operation
    # if not os.path.isfile(fileout+'.ecsv') or update:
    if True: # always update!
        # read in contoursdict
        # measurements
        contours, xc, yc=[contoursdict[x] for x in ['contours','xc','yc']]
        dictisoshape=iso.ShapeParamsdict_from_contours(contours, xc, yc)
        kwargs={'pixelsize':0.396, 'pixunit':u.Unit('arcsec')}
        tabisoshape=makeParamsuTable_from_Dict(dictisoshape, useunits=True, **kwargs)

        if towritetab:
            # construct header
            hkeys=['fileimg', 'filecontoursdict', 'isocut_rest', 'isoareallimit', 'isothreshold', 'xc', 'yc', 'masked', 'mask']
            tabheader=at.Table()
            for key in hkeys:
                if isinstance(contoursdict[key], u.quantity.Quantity):
                    tabheader[key]=[str(contoursdict[key])] 
                else:
                    tabheader[key]=[contoursdict[key]]
                if key == 'isothreshold':
                    tabheader[key]=[round(contoursdict[key],3)]

            # construct measurement
            for col in tabisoshape.colnames:
                tabisoshape.rename_column(col,'isoshape_'+col)
            tabrow=at.hstack([tabheader, tabisoshape])
            # writing
            if os.path.isfile(fileout+'.ecsv'):
                tabold=at.Table.read(fileout+'.ecsv',format='ascii.ecsv')
                tabout=at.vstack([tabold,tabrow])
                # remove duplicated old rows from table
                nr=len(tabout)-1
                hnew=at.Table(tabout[hkeys][nr])
                hold=at.Table(tabout[hkeys][:nr])
                if hnew in hold:
                    tabout.remove_rows(np.where(hold==hnew)[0])
            else: 
                tabout=tabrow
            tabout.write(fileout+'.ecsv',format='ascii.ecsv', overwrite=True)
            tabout.write(fileout+'.csv',format='ascii.csv', overwrite=True)
        if toplot:
            plotIsoMsr_fromContoursDict(dir_obj, fileplot+'.pdf', contoursdict=contoursdict, dictisoshape=dictisoshape)
    else: 
        print "skip dir_MeasureImgIso_fromContoursDict() as file exists "+fileout+'.ecsv'

    
def dir_MeasureImgIso(dir_obj, fileimg='stamp-lOIII5008_I_norm.fits', filecontoursdict=None, isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'), isoareallimit=0, towritetab=True, toplot=True, update=True):
    """
    Make iso measurements on the specified image in the directory, and write 
    the results in file dir_obj+'measureimg.csv'. 
    
    Parameters
    ------
    dir_obj: string
        path to directory
    fileimg='stamp-lOIII5008_I_norm.fits': string
        the image to measure
    isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)')
        the rest frame isophotal cut
    isoareallimit: int
        the minimum area in pix that a non centroid isophote will be considered
    towritetab=True, toplot=True, update=False

    Returns
    ------
    None

    Write Output
    ------
    dir_obj+'measureimg.csv'

    Note
    -----
    params dictionary do not have physical units, all are in pixel native units, while params table have physical units. 
    """
    # set up filenames
    filein = dir_obj+fileimg
    # fileout = dir_obj+'measureimg_iso'
    # start operation
    if True:  # always update!
        contoursdict = dir_makeContoursDict(dir_obj, fileimg=fileimg, filecontoursdict=filecontoursdict, isocut_rest=isocut_rest, isoareallimit=isoareallimit, update=update)

        dir_MeasureImgIso_fromContoursDict(dir_obj, contoursdict, towritetab=towritetab, toplot=toplot, update=update)
    else: 
        print "skip dir_MeasureImgIso() as file exist"


def plotIsoMsr(fileplot, img, threshold, contours, dictisoshape=None):
    """
    dictisoshape: 
        dictionary with keys 
        [dferetmax, theta_dferetmax, rferetmax, theta_rferetmax]

    """
    # xc, yc = 0.5*(np.array(img.shape)-1.)
    xc, yc = standards.get_img_xycenter(img)

    allcontours=polytools.find_isocontours(img,threshold)

    plt.clf()
    fig,ax=plottools.plot_img(img,vmin=0.,vmax=5.*threshold)
    if len(allcontours)>0:
        plottools.overplot_contours(ax, allcontours,color='grey',lw=2)
    if len(contours)>0:
        plottools.overplot_contours(ax, contours,color='white')
        if dictisoshape is not None:
            plottools.overplot_feretDR_frommdict(ax, dictparam=dictisoshape, xc=xc, yc=yc, color='black')
    fig.savefig(fileplot)


def plotIsoMsr_fromContoursDict(dir_obj, fileplot, contoursdict, dictisoshape=None):
    """
    dictisoshape: 
        dictionary with keys 
        [dferetmax, theta_dferetmax, rferetmax, theta_rferetmax]

    """
    contours, fileimg, isothreshold, xc, yc=[contoursdict[x] for x in ['contours', 'fileimg', 'isothreshold', 'xc', 'yc',]]
    img=fits.getdata(dir_obj+fileimg)

    plotIsoMsr(fileplot, img, isothreshold, contours, dictisoshape=dictisoshape)


def makeParamsuTable_from_Dict(paramsdict, useunits=True, pixelsize=0.396, pixunit=u.Unit('arcsec'), imgunit=1.e-15*u.Unit('erg s-1 cm-2 arcsec-2')):
    """
    transform params dictionary to table with the correct units conversion. 
    """
    # setting
    pixscale= pixelsize * (pixunit)
    scalings={'theta': u.Unit('degree'),
              'flux': imgunit*pixscale**2 ,
              'area': (pixscale**2),
              'xc': pixscale, 'yc': pixscale,
              'a': pixscale, 'b': pixscale,
              'dferetmax': pixscale, 'rferetmax': pixscale, 
              'dferetper': pixscale, 
              'theta_dferetmax': u.Unit('degree'),
              'theta_rferetmax': u.Unit('degree'),
              'theta_dferetper': u.Unit('degree')
              }

    # operation
    paramvalues=paramsdict.values()
    paramkeys=paramsdict.keys()

    if not useunits:
        paramlist=[[value] for value in paramvalues]
        tabout=at.Table(paramlist,names=paramkeys)
        return tabout
    else: 
        # assign units
        paramsu = [None] * len(paramvalues) # paramvalues with units
        for i, key in enumerate(paramkeys):
            if key in scalings:
                paramsu[i]=paramvalues[i]*scalings[key]
            else:
                paramsu[i]=paramvalues[i]

        # assemble table
        tabout=at.Table()
        for i, key in enumerate(paramkeys):
            if isinstance(paramsu[i],u.quantity.Quantity):
                tabout[key]=[paramsu[i].value]
                tabout[key].unit=paramsu[i].unit
            else: 
                tabout[key]=[paramsu[i]]
        return tabout



# going to be deprecated

def measureparams(img, method='moment', sigma=0., isothreshold=0., isoareallimit=0., pixelsize=0.396, imgunit=1.e-15*u.Unit('erg s-1 cm-2 arcsec-2'), pixunit=u.Unit('arcsec'), useunits=True,  tooffset=False): 
    """
    Infer the ellipse params:
        ('flux' or 'area','xc','yc','a','b','theta','eccen')

    from image using specified methods ('moment' or 'gaussfit')
    with units and pixsize handeling (default useunits=True).  

    Its an option to offset image (minus median) by setting 
    tooffset=True. 

    SDSS pixel size (0.396 arcsec) is assuemd by default. 

    Parameters
    ------

    img: np 2d array
        line intensity map (with or without unit)
    method='moment': string
        specify what method to use
            moment: 2d moments
            gaussfit: fit image to 2d gaussian
            isofit: fit ellipse to iso contour
            isomom: 2d moments of iso contour
    sigma=0.
        image noise level for gaussfit
    isothreshold=0.
        isophotal isothreshold for iso
    pixelsize=0.396
        the image pixel size (without unit)
    imgunit=u.Unit('erg s-1 cm-2 arcsec-2')
        image unit
    pixunit=u.Unit('arcsec')
        pixel size unit
    useunits=True (bool)
        if true, output will have units as  
        specified. otherwise, no units, and 
        distances are measured in pixels. 

    RETURN: 
        tabout: astropy table of 1 row with columns: 
             flux   : total integrated flux (if method = moment or gaussfit)
                      (default unit erg / (cm2 s)) 
             area   : net enclosed area of iso contour (if method = iso)
                      (default unit arcsec2) 
             xc  : x axis centroid relative to image center
                      (default unit arcsec)
             yc  : y axis centroid relative to image center
                      (default unit arcsec) 
             a : length of the major axis (sigma)
                      (default unit arcsec)
             b : length of the minor axis (sigma)
                      (default unit arcsec)
             theta  : orientation of the major axis, y of x
                      (unit degree)
             eccen : eccentricity of the ellipse 
                            sqrt(1-(minor/major)^2)

    DESCRIPTION: 
        see moment and gaussfit modules
    """
    #==== init
    pixscale= pixelsize * (pixunit)
    scalings={'theta': u.Unit('degree'),
          'flux': imgunit*pixscale**2 ,
          'area': (pixscale**2),
          'xc': pixscale, 'yc': pixscale,
          'a': pixscale, 'b': pixscale,
          'dferetmax': pixscale, 'rferetmax': pixscale, 
          'dferetper': pixscale, 
          'theta_dferetmax': u.Unit('degree'),
          'theta_rferetmax': u.Unit('degree'),
          'theta_dferetper': u.Unit('degree')
          }


    #==== protect the original params from being altered
    from copy import deepcopy
    img=deepcopy(img)
    if tooffset: img=img-np.median(img)

    if useunits:
        if hasattr(img, 'unit'):
            if img.unit != imgunit: raise ValueError('unit inconsistent')
    else: 
        pixelsize,pixunit,imgunit=1.,1.,1.
    # if img has a unit get rid of it, for the convenience of calculation
    if hasattr(img, 'unit'): img=img.value

    #===== calculating the params
    if method=='moment':
        cols=['flux','xc','yc','a','b','theta']
        params=moment.mom_ellipseparams(img)
    elif method == 'gaussfit':
        cols=['flux','xc','yc','a','b','theta']
        if sigma==0.: raise NameError('please set sigma')
        params=gaussfit.gaussfit_ellipseparams(img,sigma)
    elif method == 'isofit':
        cols=['area','xc','yc','a','b','theta']
        if isothreshold==0.: raise NameError('please set isothreshold')
        params=iso.iso_FitEllipseParams(img,isothreshold)
    elif method == 'isomom':
        cols=['area','xc','yc','a','b','theta']
        if isothreshold==0.: raise NameError('please set isothreshold')
        params=iso.iso_MomEllipseParams(img, isothreshold)
    elif method == 'isoshape':
        if isothreshold==0.: raise NameError('please set isothreshold')
        if isoareallimit==0.: raise NameError('please set isoareallimit')
        paramsdict = iso.iso_ShapeParams(img, isothreshold, isoareallimit)
        params = paramsdict.values()
        cols = paramsdict.keys()
        # ['area', 'dferetmax', 'theta_dferetmax', 'rferetmax', 'theta_rferetmax', 'dferetper', 'theta_dferetper', 'aspectr']
    else: 
        raise NameError('method not recognized')

    
    #===== assemble the output
    if not useunits:
        paramlist=[[params[i]] for i in range(len(params))]
        tabout=at.Table(paramlist,names=cols)
        return tabout
    else: 
        # assign units
        paramsu = [None] * len(params)
        for i, var in enumerate(cols):
            if var in scalings:
                paramsu[i]=params[i]*scalings[var]
            else:
                paramsu[i]=params[i]

        # assemble table
        tabout=at.Table()
        for i, var in enumerate(cols):
            if isinstance(paramsu[i],u.quantity.Quantity):
                tabout[var]=[paramsu[i].value]
                tabout[var].unit=paramsu[i].unit
            else: 
                tabout[var]=[paramsu[i]]
        return tabout
    


