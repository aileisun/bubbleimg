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

from .. import standards


def dir_doIsos(dir_obj, isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'),  isoareallimit=10, contrastr=0.1, update=True, toplot=True): 
    """
    do the following operations: 
        find/store isocut contour on denoised line map as linblob
        measure isoshape on the lineblob
        find/store isocut/contrastr contour on continuum map as galmask
        subtract galmask from lineblob as lineblobmasked
        measure isoshape on the lineblobmasked
    """
    fileimg_line = 'stamp-lOIII5008_I_norm_denoised.fits'
    fileimg_line_psfresid = 'stamp-lOIII5008_I_norm_psfresidual_denoised.fits'
    fileimg_gal = 'stamp-conti-onOIIIscale_I_norm.fits'
    fileimg_psfmodel = 'stamp-lOIII5008_I_norm_psfmodel.fits'
    
    prefix_cdict='contoursdict_'
    fcdict_linblob = prefix_cdict+'linblob.pkl'
    fcdict_galmask = prefix_cdict+'galmask.pkl'
    fcdict_galctrmask = prefix_cdict+'galctrmask.pkl'
    fcdict_linblob_galm = prefix_cdict+'linblob_galmasked.pkl'
    fcdict_psfmask = prefix_cdict+'psfmask.pkl'
    fcdict_linblob_galm_psfm = prefix_cdict+'linblob_galmasked_psfmasked.pkl'
    fcdict_linblob_psfresid =  prefix_cdict+'linblob_psfresid.pkl'
    fcdict_linblob_psfresid_galm =  prefix_cdict+'linblob_psfresid_galmasked.pkl'

    # # start temporary codes
    # cdict_galmask = read_pickle(dir_obj+fcdict_galmask)
    # cdict_linblob = read_pickle(dir_obj+fcdict_linblob)

    # # end temporary codes

    # get contour of line blob
    cdict_linblob=dir_makeContoursDict(dir_obj, fileimg=fileimg_line, filecontoursdict=fcdict_linblob,isocut_rest=isocut_rest, isoareallimit=isoareallimit, update=update)

    # get contour of galaxy mask
    cdict_galmask=dir_makeContoursDict(dir_obj, fileimg=fileimg_gal, filecontoursdict=fcdict_galmask, isocut_rest=isocut_rest/contrastr, isoareallimit=0, update=update)

    cdict_galmaskctr = dir_makeContoursDict(dir_obj, fileimg=fileimg_gal, filecontoursdict=fcdict_galctrmask, isocut_rest=isocut_rest/contrastr, isoareallimit=64*64, update=update)


    # measure isoshape of lineblob
    dir_MeasureImgIso_fromContoursDict(dir_obj, cdict_linblob, towritetab=True, toplot=toplot, update=update)

    # plot isoshape of galmask
    dir_MeasureImgIso_fromContoursDict(dir_obj, cdict_galmask, towritetab=False, toplot=toplot, update=update)

    # mask galmask from linblob -> linblob_galm
    if not os.path.isfile(dir_obj+fcdict_linblob_galm) or update:
        cdict_linblob_galm = shapelytools.diffcontoursdict(cdict_linblob, cdict_galmask, fcdict_linblob_galm)
        write_pickle(cdict_linblob_galm, dir_obj+fcdict_linblob_galm)

    # measure isoshape from lineblob_galm
    dir_MeasureImgIso_fromContoursDict(dir_obj, cdict_linblob_galm, towritetab=True, toplot=toplot, update=update)

    # measure iso for psfmodel to make psfmask -- optional
    kwargs = {'fileimg': fileimg_psfmodel,
              'filecontoursdict':fcdict_psfmask, 
              'isocut_rest': isocut_rest,
              'isoareallimit': 0,
              'update': update}
    dir_MeasureImgIso(dir_obj, **kwargs)

    # measure iso for psf residual 
    kwargs = {'fileimg': fileimg_line_psfresid,
              'filecontoursdict':fcdict_linblob_psfresid, 
              'isocut_rest': isocut_rest,
              'isoareallimit': isoareallimit,
              'update': update}
    dir_MeasureImgIso(dir_obj, **kwargs)


    # mask galmask from linblob_psfresid -> linblob_psfresid_galm
    cdict_linblob_psfresid = read_pickle(dir_obj+fcdict_linblob_psfresid)
    if not os.path.isfile(dir_obj+fcdict_linblob_psfresid_galm) or update:
        cdict_linblob_psfresid_galm = shapelytools.diffcontoursdict(cdict_linblob_psfresid, cdict_galmask, fcdict_linblob_psfresid_galm)
        write_pickle(cdict_linblob_psfresid_galm, dir_obj+fcdict_linblob_psfresid_galm)

    # measure iso for psf residual galm
    cdict_linblob_psfresid_galm = read_pickle(dir_obj+fcdict_linblob_psfresid_galm)
    dir_MeasureImgIso_fromContoursDict(dir_obj, cdict_linblob_psfresid_galm, towritetab=True, toplot=toplot, update=update)



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
            tabout.write(fileout+'.ecsv',format='ascii.ecsv')
            tabout.write(fileout+'.csv',format='ascii.csv')
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
    scalings={'theta':u.Unit('degree'),
              'flux': imgunit*pixscale**2 ,
              'area':(pixscale**2),
              'xc':pixscale,'yc':pixscale,
              'a':pixscale,'b':pixscale,
              'dferetmax':pixscale,'rferetmax':pixscale,
              'theta_dferetmax':u.Unit('degree'),
              'theta_rferetmax':u.Unit('degree'),
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
    scalings={'theta':u.Unit('degree'),
              'flux': imgunit*pixscale**2 ,
              'area':(pixscale**2),
              'xc':pixscale,'yc':pixscale,
              'a':pixscale,'b':pixscale,
              'dferetmax':pixscale,'rferetmin':pixscale,
              'theta_dferetmax':u.Unit('degree'),
              'theta_rferetmin':u.Unit('degree'),
              # 'feretmax':pixscale,'feretmin':pixscale,
              # 'feret90':pixscale,
              # 'theta_feretmax':u.Unit('degree'),
              # 'theta_feretmin':u.Unit('degree'),
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
        # cols=['area','feretmax','theta_feretmax','feretmin',
        #       'theta_feretmin','feret90','aspectratio']
        cols=['area', 'dferetmax', 'theta_dferetmax', 'rferetmax', 'theta_rferetmax']
        if isothreshold==0.: raise NameError('please set isothreshold')
        if isoareallimit==0.: raise NameError('please set isoareallimit')
        params=iso.iso_ShapeParams(img, isothreshold, isoareallimit)
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
    



# Deprecated functions:
# def plotIsoMsr(fileplot, img, threshold, tabisoshape):
#     ... to be revised ...
#     # plot iso measurements
#     maincontour, holecontours=iso.find_centercontour(img, threshold=threshold)
#     contours=iso.find_contours(img, threshold=threshold)
#     plt.clf()
#     fig,ax=plottools.plot_img(img,vmin=0.,vmax=threshold*5.)

#     if len(contours)>0:
#         plottools.overplot_contours(ax, contours,color='grey',lw=2)
#     if maincontour is not None:
#         plottools.overplot_contours(ax, [maincontour]+holecontours,color='white')
#         plottools.overplot_ellipse_fromtab(ax,img,tabisofit,color='red',useunits=True)
#         plottools.overplot_ellipse_fromtab(ax,img,tabisomom,color='green',useunits=True)
#         plottools.overplot_feret_fromtab(ax,img,tabisoshape,tabisomom,color='blue',useunits=True)
#     fig.savefig(fileplot)

# def dir_MakeContourPlot_fits(dir_obj,filename='stamp-lOIII5008_I_norm.fits',isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'),toplot=True):
#     """ make iso contour plot """

#     fileimg_base=os.path.splitext(filename)[0]
#     filein=dir_obj+filename

#     # read in 
#     img=fits.getdata(filein)
#     header=fits.getheader(filein)
    
#     # calculate iso threshold
#     xid=Table.read(dir_obj+'xid.csv',format='ascii.csv',comment='#')
#     z=xid['z'][0]
#     imgunit=u.Unit(header['BUNIT'])
#     isocut_obs=isocut_rest*(1.+z)**-4
#     threshold=(isocut_obs/imgunit).to(u.dimensionless_unscaled).value


#     # plot iso measurements
#     fileplot=dir_obj+fileimg_base+'_isocontour.pdf'
#     maincontour, holecontours =iso.find_centercontour(img,threshold=threshold)
#     plt.clf()
#     fig,ax=plottools.plot_img(img,vmin=0.,vmax=threshold*5.)
#     if maincontour is not None:
#         plottools.overplot_contours(ax, [maincontour]+holecontours,color='white')
#     fig.savefig(fileplot)

#     # plot iso measurements
#     fileplot=dir_obj+fileimg_base+'.pdf'
#     maincontour, holecontours=iso.find_centercontour(img,threshold=threshold)
#     plt.clf()
#     fig,ax=plottools.plot_img(img,vmin=0.,vmax=threshold*5.)
#     fig.savefig(fileplot)

# def dir_MeasureImgIso(dir_obj,filename='stamp-lOIII5008_I.fits',isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'),toplot=True,update=False):
#     """
#     Make iso measurements on the specified image in the directory, and write 
#     the results in file dir_obj+'measureimg.csv'. 
    
#     Parameters
#     ------
#     dir_obj: string
#         path to directory
#     filename='stamp-lOIII5008_I.fits': string
#         the image to measure

#     Returns
#     ------
#     tabout: astropy Table of measurements

#     Write Output
#     ------
#     dir_obj+'measureimg.csv'
#     """
#     # set up
#     fileimg_base=os.path.splitext(filename)[0]
#     filein=dir_obj+filename
#     fileout=dir_obj+'measureimg'+'_'+fileimg_base
#     fileplot=dir_obj+fileimg_base+'_iso.pdf'

#     if not os.path.isfile(fileout+'.ecsv') or update:
#         # read in 
#         img=fits.getdata(filein)
#         header=fits.getheader(filein)
        
#         # calculate iso threshold
#         xid=Table.read(dir_obj+'xid.csv',format='ascii.csv',comment='#')
#         z=xid['z'][0]
#         imgunit=u.Unit(header['BUNIT'])
#         isocut_obs=isocut_rest*(1.+z)**-4
#         threshold=(isocut_obs/imgunit).to(u.dimensionless_unscaled).value


#         # calculation
#         kwargs={'threshold':threshold,'pixelsize':0.396,'imgunit':imgunit,'pixunit':u.Unit('arcsec'),'useunits':True,'tooffset':False}

#         tabisofit=measureparams(img,method='isofit',**kwargs)
#         tabisomom=measureparams(img,method='isomom',**kwargs)
#         tabisoshape=measureparams(img,method='isoshape',**kwargs)

#         if toplot:
#             plotmsr(fileplot, img, threshold, tabisofit, tabisomom, tabisoshape)

#         # assemble table
#         for col in tabisofit.colnames:
#             tabisofit.rename_column(col,'isofit_'+col)
#         for col in tabisomom.colnames:
#             tabisomom.rename_column(col,'isomom_'+col)
#         for col in tabisoshape.colnames:
#             tabisoshape.rename_column(col,'isoshape_'+col)

#         tabheader=at.Table([[filename],[z],[str(isocut_rest)],[threshold]],names=['filename','z','isocut_rest','threshold'])
#         tabout=at.hstack([tabheader,tabisofit,tabisomom,tabisoshape])
#         # writing
#         tabout.write(fileout+'.ecsv',format='ascii.ecsv')
#         tabout.write(fileout+'.csv',format='ascii.csv')
#         # tabout.write(fileout+'.fits',format='fits')
#     else: 
#         print "skip dir_MeasureImgIso() as file exist"

# def plotmsr(fileplot, img, threshold, tabisofit, tabisomom, tabisoshape):
#     # plot iso measurements
#     maincontour, holecontours=iso.find_centercontour(img, threshold=threshold)
#     contours=iso.find_contours(img, threshold=threshold)
#     plt.clf()
#     fig,ax=plottools.plot_img(img,vmin=0.,vmax=threshold*5.)

#     if len(contours)>0:
#         plottools.overplot_contours(ax, contours,color='grey',lw=2)
#     if maincontour is not None:
#         plottools.overplot_contours(ax, [maincontour]+holecontours,color='white')
#         plottools.overplot_ellipse_fromtab(ax,img,tabisofit,color='red',useunits=True)
#         plottools.overplot_ellipse_fromtab(ax,img,tabisomom,color='green',useunits=True)
#         plottools.overplot_feret_fromtab(ax,img,tabisoshape,tabisomom,color='blue',useunits=True)
#     fig.savefig(fileplot)



    # # plot the galmask
    # fileout=dir_obj+'measureimg_iso'
    # fileplot=fileout+'_'+os.path.splitext(fcdict_galmask)[0]+'.pdf'
    # if not os.path.isfile(fileplot) or update:
    #     if toplot:
    #         plotIsoMsr_fromContoursDict(dir_obj, fileplot, cdict_galmask, dictisoshape=None)


    # # mask psfmask from linblob_galm
    # cdict_psfmask=read_pickle(dir_obj+fcdict_psfmask)

    # if not os.path.isfile(dir_obj+fcdict_linblob_galm_psfm) or update:
    #     cdict_linblob_galm_psfm=shapelytools.diffcontoursdict(cdict_linblob_galm, cdict_psfmask, fcdict_linblob_galm_psfm)
    #     write_pickle(cdict_linblob_galm_psfm, dir_obj+fcdict_linblob_galm_psfm)
    # else:
    #     cdict_linblob_galm_psfm = read_pickle(dir_obj+fcdict_linblob_galm_psfm)

    # # measure iso for psf masked galaxy masked line blob
    # dir_MeasureImgIso_fromContoursDict(dir_obj, cdict_linblob_galm_psfm, towritetab=True, toplot=True, update=update)
