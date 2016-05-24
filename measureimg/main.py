# main.py
# ALS 2016/04/27

"""
Interface with outer funcitons to handle all image measurements (moment, gaussfit) with specifics like units and pixel sizes. 

"""
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.table import Table, hstack
from astropy.io import fits
import os

import moment
import gaussfit
import iso
import plottools

# import sys
# sys.path.append('../denoiseimg')
# import noiselevel
# sigma=noiselevel.load_noiselevel(dir_obj,filename,update=False)



def dir_MakeContourPlot_fits(dir_obj,filename='stamp-lOIII5008_I_norm.fits',isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'),toplot=True):
    """ make contour plot """

    filename_base=os.path.splitext(filename)[0]
    filein=dir_obj+filename

    # read in 
    img=fits.getdata(filein)
    header=fits.getheader(filein)
    
    # calculate iso threshold
    xid=Table.read(dir_obj+'xid.csv',format='ascii.csv',comment='#')
    z=xid['z'][0]
    imgunit=u.Unit(header['BUNIT'])
    isocut_obs=isocut_rest*(1.+z)**-4
    threshold=(isocut_obs/imgunit).to(u.dimensionless_unscaled).value


    # plot iso measurements
    fileplot=dir_obj+filename_base+'_isocontour.pdf'
    maincontour, holecontours=iso.find_centercontour(img,threshold=threshold)
    plt.clf()
    fig,ax=plottools.plot_img(img,vmin=0.,vmax=threshold*5.)
    if maincontour is not None:
        plottools.overplot_contours(ax, [maincontour]+holecontours,color='white')
    fig.savefig(fileplot)

    # plot iso measurements
    fileplot=dir_obj+filename_base+'.pdf'
    maincontour, holecontours=iso.find_centercontour(img,threshold=threshold)
    plt.clf()
    fig,ax=plottools.plot_img(img,vmin=0.,vmax=threshold*5.)
    fig.savefig(fileplot)



def dir_MeasureImgIso_fits(dir_obj,filename='stamp-lOIII5008_I.fits',isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'),toplot=True):
    """
    Make iso measurements on the specified image in the directory, and write 
    the results in file dir_obj+'measureimg.csv'. 
    
    Parameters
    ------
    dir_obj: string
        path to directory
    filename='stamp-lOIII5008_I.fits': string
        the image to measure

    Returns
    ------
    tabout: astropy Table of measurements

    Write Output
    ------
    dir_obj+'measureimg.csv'
    """
    # set up
    filename_base=os.path.splitext(filename)[0]
    filein=dir_obj+filename
    fileout=dir_obj+'measureimg'+'_'+filename_base
    fileplot=dir_obj+filename_base+'_iso.pdf'

    # read in 
    img=fits.getdata(filein)
    header=fits.getheader(filein)
    
    # calculate iso threshold
    xid=Table.read(dir_obj+'xid.csv',format='ascii.csv',comment='#')
    z=xid['z'][0]
    imgunit=u.Unit(header['BUNIT'])
    isocut_obs=isocut_rest*(1.+z)**-4
    threshold=(isocut_obs/imgunit).to(u.dimensionless_unscaled).value


    # calculation
    kwargs={'threshold':threshold,'pixelsize':0.396,'imgunit':imgunit,'pixunit':u.Unit('arcsec'),'useunits':True,'tooffset':False}

    tabisofit=measureparams(img,method='isofit',**kwargs)
    tabisomom=measureparams(img,method='isomom',**kwargs)
    tabisoshape=measureparams(img,method='isoshape',**kwargs)

    if toplot:
        # plot iso measurements
        maincontour, holecontours=iso.find_centercontour(img,threshold=threshold)
        plt.clf()
        fig,ax=plottools.plot_img(img,vmin=0.,vmax=threshold*5.)
        if maincontour is not None:
            plottools.overplot_contours(ax, [maincontour]+holecontours,color='white')
            plottools.overplot_ellipse_fromtab(ax,img,tabisofit,color='red',useunits=True)
            plottools.overplot_ellipse_fromtab(ax,img,tabisomom,color='green',useunits=True)
            plottools.overplot_feret_fromtab(ax,img,tabisoshape,tabisomom,color='blue',useunits=True)
        fig.savefig(fileplot)


    # assemble table
    for col in tabisofit.colnames:
        tabisofit.rename_column(col,'isofit_'+col)
    for col in tabisomom.colnames:
        tabisomom.rename_column(col,'isomom_'+col)
    for col in tabisoshape.colnames:
        tabisoshape.rename_column(col,'isoshape_'+col)

    tabheader=Table([[filename],[z],[str(isocut_rest)],[threshold]],names=['filename','z','isocut_rest','threshold'])
    # tabheader=Table([[xid['ra'][0]],[xid['dec'][0]],[z],[filename],[str(isocut_rest)],[threshold]],names=['RA','DEC','z','filename','isocut_rest','threshold'])
    tabout=hstack([tabheader,tabisofit,tabisomom,tabisoshape])
    # writing
    tabout.write(fileout+'.ecsv',format='ascii.ecsv')
    tabout.write(fileout+'.csv',format='ascii.csv')
    # tabout.write(fileout+'.fits',format='fits')



def measureparams(img,method='moment',sigma=0.,threshold=0.,pixelsize=0.396,imgunit=u.Unit('erg s-1 cm-2 arcsec-2'),pixunit=u.Unit('arcsec'),useunits=True, tooffset=False): 
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
    threshold=0.
        isophotal threshold for iso
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
    scalings={'flux': imgunit*pixscale**2 ,
              'area':(pixscale**2),
              'xc':pixscale,'yc':pixscale,
              'a':pixscale,'b':pixscale,
              'feretmax':pixscale,'feretmin':pixscale,
              'feret90':pixscale,
              'theta':u.Unit('degree'),
              'theta_feretmax':u.Unit('degree'),
              'theta_feretmin':u.Unit('degree'),
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
        if threshold==0.: raise NameError('please set threshold')
        params=iso.iso_FitEllipseParams(img,threshold)
    elif method == 'isomom':
        cols=['area','xc','yc','a','b','theta']
        if threshold==0.: raise NameError('please set threshold')
        params=iso.iso_MomEllipseParams(img,threshold)
    elif method == 'isoshape':
        cols=['area','feretmax','theta_feretmax','feretmin',
              'theta_feretmin','feret90','aspectratio']
        if threshold==0.: raise NameError('please set threshold')
        params=iso.iso_ShapeParams(img,threshold)
    else: 
        raise NameError('method not recognized')

    
    #===== assemble the output
    if not useunits:
        paramlist=[[params[i]] for i in range(len(params))]
        tabout=Table(paramlist,names=cols)
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
        tabout=Table()
        for i, var in enumerate(cols):
            if isinstance(paramsu[i],u.quantity.Quantity):
                tabout[var]=[paramsu[i].value]
                tabout[var].unit=paramsu[i].unit
            else: 
                tabout[var]=[paramsu[i]]
        return tabout
    


    # #===== adding in the units
    # if method in ['moment','gaussfit']:
    #     flux=flux*pixelsize**2 * (imgunit*pixunit**2)
    # elif method in ['isofit','isomom']:
    #     area=area*pixelsize**2
    # xc = xc*pixelsize * (pixunit)
    # yc = yc*pixelsize * (pixunit)
    # a=a*pixelsize * (pixunit)
    # b=b*pixelsize * (pixunit)
    # theta=theta
    # eccen = np.sqrt(1.-b/a)



    # # # assign the variables
    # for i, var in enumerate(cols):
    #     exec(var + " = params["+str(i)+"]")

    # if method in ['moment','gaussfit','isofit','isomom']:
    #     cols=cols+['aspectratio']
    #     aspectratio=b/a
