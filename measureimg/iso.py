# iso.py
# ALS 2015/04/22

"""
Make and return isophotal measurements
"""

import numpy as np
import astropy.units as u
from astropy.table import Table

import polytools

import collections

def ShapeParamsdict_from_contours(contours, xc, yc):
    """
    return a dictionary of key:value pairs of ShapeParams measured from contours
    """
    cols=['area', 'dferetmax', 'theta_dferetmax', 'rferetmax', 'theta_rferetmax', 'dferetper', 'theta_dferetper', 'aspectr']
    if len(contours)>0:
        area = polytools.NetPolygonsArea(contours)
        dferetmax, theta_dferetmax = polytools.FeretD_max(contours)
        rferetmax, theta_rferetmax = polytools.FeretR_max(contours, xc, yc)
        dferetper, theta_dferetper = polytools.FeretD_max90(contours)
        aspectr = polytools.FeretAspectRatio(contours)
        values=[area, dferetmax, theta_dferetmax, rferetmax, theta_rferetmax, dferetper, theta_dferetper, aspectr]
    else: 
        values=[0., 0., np.nan, 0., np.nan, 0., np.nan, np.nan]

    paramsdict=collections.OrderedDict(zip(cols, values))
    return paramsdict

# def iso_ShapeParams(img, threshold, areallimit):
#     """
#     Make shape measurements to the iso contours. The contour at the center 
#     and the contours larger than an area of areallimit are considered. 

#     Parameters
#     ------
#     img (2d np array)
#     threshold (float): threshold of the isophote

#     Returns
#     ------
#     params: list
#         area: net area in the centroid "isophote" (minus the holes). 
#         dferetmax: maximum feret's diameter
#         theta_dferetmax: orientation of the dferetmax in degree y of x

#         rferetmax: maximum radius measured from center
#         theta_rferetmax: orientation of the rferetmax in degree y of x

#     DESCRIPTION
#     ------
#     see polytools
#     """
#     xc, yc = 0.5*(np.array(img.shape)-1.)

#     # find contour
#     contours=polytools.find_realisocontours(img, threshold, areallimit)

#     # calculate
#     if len(contours)>0:
#         area=polytools.NetPolygonsArea(contours)
#         dferetmax, theta_dferetmax=polytools.FeretD_max(contours)
#         rferetmax, theta_rferetmax=polytools.FeretR_max(contours, xc, yc)

#         return [area, dferetmax, theta_dferetmax, rferetmax, theta_rferetmax]
#     else: 
#         return [0, 0, 0, 0, 0]




# def iso_MomEllipseParams(img, threshold):
#     """
#     Get the moments (represented in ellipse params) of the region enclosed by 
#     the center isophote at threshold. Its an estimate as the region is 
#     pixelized. Holes inside the region is not taken into accoutn (yet). 

#     Parameters
#     ------
#     img (2d np array)
#     threshold (float): threshold of the isophote

#     Returns
#     ------
#     params: list
#         area: area enclosed in the isophote, holes are not subtracted
#         xc: x centroid
#         yc: y centroid
#         a (length in pixels): sigma_major 
#         b (length in pixels): sigma_minor
#         theta (orientation in degree): y of x

#     DESCRIPTION
#     ------
#     uses polytools.mom_ellipseparams_poly(poly)
#     """
#     maincontour, holecontours = polytools.find_centercontour(img,threshold)

#     if maincontour is not None:
#         params=polytools.mom_ellipseparams_poly(maincontour)
#         return params
#     else:
#         return [0.,0.,0.,0.,0.,0.]


# def iso_FitEllipseParams(img, threshold):
#     """
#     Get params of the best fit ellipse to the center isophote at threshold. 

#     Parameters
#     ------
#     img (2d np array)
#     threshold (float): threshold of the isophote

#     Returns
#     ------
#     params: list
#         area: area of the ellipse, pi * a * b
#         xc: x centroid
#         yc: y centroid
#         a (length in pixels): semi-major axis
#         b (length in pixels): semi-minor axis
#         theta (orientation in degree): y of x

#     DESCRIPTION
#     ------
#     see fitEllipsetoContour() and polytools.find_centercontour()
#     """
#     maincontour, holecontours = polytools.find_centercontour(img,threshold)

#     if maincontour is not None:
#         params=list(fitEllipsetoContour(maincontour)) # yc, xc, a, b, theta
#         ellipsearea=np.pi*params[2]*params[3]
#         return [ellipsearea]+params
#     else:
#         return [0.,0.,0.,0.,0.,0.]



def fitEllipsetoContour(poly):
    """
    Fit an ellipse to a polygon

    Parameters
    ----------
    poly : (N*2) ndarray of double
        list of corner coordinates (x,y) of the polygon. N is the number of 
        corners. 

    Returns
    -------
    params : list (xc, yc, a, b, theta)
            xc: centroid x position
            yc: centroid y position
            a:  semi-major axis
            b:  semi-minor axis  (all in pix)
            theta: orientation of semi-major axis y of x in degrees

    Note
    -------
    uses skimage.measure.EllipseModel. Note that it has a wierd convention of 
    (row, col), where row is y and col is x. So the order is reversed. 
    EllipseModel() can also return predicted model of ellipse (as a polygon), 
    and the residual. But these functions are not used. 
    # modelcontour=model.predict_xy(np.linspace(0,2.*np.pi))

    Reference
    -------
    http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.EllipseModel

    """
    from skimage.measure import EllipseModel

    model = EllipseModel()
    success=model.estimate(poly)

    if success:
        yc, xc, a, b, theta=model.params # xc, yc, a, b, theta <- in radiand
        # convert theta to degree
        theta=np.degrees(theta)
        theta=90.-theta
        # make sure major is larger than minor and theta is for major
        if a<b:
            [a,b]=[b,a]
            theta=theta-90.
        if theta<0.: theta=theta+180.
        params=(xc, yc, a, b, theta)
        return params
    else:
        return (0, 0, 0, 0, 0)



# def find_centercontour_old(img,threshold):
    # """
    # Find the high contour at the of the image, and all the low contours 
    # enclosed in it. 

    # Parameters
    # ----------
    # img: 2d np array
    # threshold: float

    # Returns
    # -------
    # maincontour: (N*2) ndarray of double
    #     the polygon of the center high contour. points are counter-clockwise.
    #     If failed, return []. 
    # holecontours:
    #     a list of the polygons of the holes (low contour) enclosed in the 
    #     center contours. Points are clockwise. 
    # """
    # from skimage import skimage.measure

    # contours=find_contours(img,threshold)
    
    # # find center positive contour which enclosed the image center
    # ccontours=[] # for storing contour(s) that are right at the center
    # xc,yc=0.5*(np.array(img.shape)-1.)
    # for i, contour in enumerate(contours):
    #     if skimage.measure.points_in_poly([[xc,yc]],contour) and polytools.SignedPolygonArea(contour)>0. :
    #         ccontours=ccontours+[contour]

    # if len(ccontours)==1: # one unique contour found
    #     maincontour = ccontours[0]

    # elif len(ccontours)>=2: # pick the smallest contour
    #     nc=len(ccontours)
    #     print "get "+str(nc)+" high center contours, pick smallest one"
    #     polyareas=np.array([polytools.SignedPolygonArea(ccontours[i]) for i in range(nc)])
    #     ic=np.argmin(polyareas)
    #     maincontour=ccontours[ic]

    #     # sanity check
    #     for i in range(nc):
    #         if i != ic:
    #             a=np.unique(skimage.measure.points_in_poly(maincontour,ccontours[i]))
    #             if len(a)==1:
    #                 if a[0]==False:
    #                     raise ValueError("the maincontour is not the smallest")
    #             else: 
    #                 raise ValueError("the two contours intersect")
    # elif len(ccontours)==0:
    #     # failed. nothing found
    #     maincontour = None
    # else: 
    #     raise ValueError('something went wrong with the countours: '+"number of ccontours:"+str(len(ccontours)))

    # # find holes in the center contour
    # if maincontour is not None:
    #     holecontours=[]
    #     for i, contour in enumerate(contours):
    #         if not np.array_equal(contour,maincontour) and np.any(skimage.measure.points_in_poly(contour,maincontour)):
    #             holecontours=holecontours+[contour]
    # else: holecontours=[]

    # return maincontour, holecontours





# def iso_ShapeParams_old(img,threshold):

#     """
#     Make shape measurements to the center isophote at threshold. 

#     Parameters
#     ------
#     img (2d np array)
#     threshold (float): threshold of the isophote

#     Returns
#     ------
#     params: list
#         area: net area in the centroid "isophote" (minus the holes). 
#         feretmax: maximum feret's diameter
#         theta_feretmax: orientation of the feretmax in degree y of x
#         feretmin: minimum feret's diameter
#         theta_feretmin: orientation of the feretmin in degree y of x
#         feret90: feret's diameter perpendicular to the maximum
#         feretaspectratio: feret90/feretmax

#     DESCRIPTION
#     ------
#     see polytools
#     """
#     # find contour
#     maincontour, holecontours = polytools.find_centercontour(img,threshold)

#     # calculate
#     if maincontour is not None:
#         area=polytools.NetPolygonsArea([maincontour]+holecontours)
#         feretmax,theta_feretmax=polytools.FeretD_max(maincontour)
#         feretmin,theta_feretmin=polytools.FeretD_min(maincontour)
#         feret90,theta_feretmax90=polytools.FeretD_max90(maincontour)
#         feretaspectratio=polytools.FeretAspectRatio(maincontour)

#         return [area,feretmax,theta_feretmax,feretmin,theta_feretmin,feret90,feretaspectratio]
#     else: 
#         return [0,0,0,0,0,0,0]



# def sandbox():
#     poly=maincontour
#     plt.close('all')
#     plt.plot(poly[0:21,0],poly[0:21,1])
#     plt.plot(poly[20:41,0],poly[20:41,1])
#     plt.plot(poly[40:61,0],poly[40:61,1])
#     plt.plot(poly[60:81,0],poly[60:81,1])
#     plt.plot(poly[80:101,0],poly[80:101,1])
#     plt.plot(poly[100:121,0],poly[100:121,1])
#     plt.show(block=False)


# def sandbox_plotting(img,threshold):
#     # filein='coif1/threshold_1sigma/img_denoised_soft.npy'
#     # img=np.load(filein)

#     maincontour, holecontours = polytools.find_centercontour(img,threshold)

#     params=fitEllipsetoContour(maincontour)

#     # import plottools
#     plt.close('all')
#     fig,ax=plottools.plot_img(img)
#     plottools.overplot_contour(ax,maincontour,color='red')
#     # plottools.overplot_contour(ax,modelcontour)
#     plottools.overplot_ellipse(ax,params,color='green')
#     fig.show()


# examples
# poly = [(2.0, 1.0), (4.0, 5.0), (7.0, 8.0)]

# def plottins():
#     # plotting
#     n=57
#     plt.close()
#     plt.clf()
#     fig=plt.figure()
#     ax=fig.add_subplot(1,1,1)
#     ax.imshow(img,interpolation='nearest',origin='lower left',aspect='equal',cmap='viridis')
#     ax.plot(maincontour[:, 1], maincontour[:, 0], linewidth=2)

#     # for i, contour in enumerate(contours):
#     #     ax.plot(contour[:, 1], contour[:, 0], linewidth=2)
#     plt.xlim(0,n)
#     plt.ylim(0,n)
#     plt.show(block=False)







# def measureiso(img,isocut,pixelsize=0.396,imgunit=u.Unit('erg s-1 cm-2 arcsec-2'),pixunit=u.Unit('arcsec'),useunits=True):
#     """
#     PURPOSE: Given isocut, measure the nebular size, area, and flux. 
#              SDSS pixel size is assuemd by default. 

#     PARAMETERS: 
#             img (2d array)     :  line intensity map [erg/s/cm2/arcsec2]
#             isocut (float):  isophoto cut [erg/s/cm2/arcsec2]
#             pixelsize=0.396    :  the image pixel size [arcsec]
#             useunits=True (bool) : if true, output will have units as  
#                                     specified. otherwise, no units, and 
#                                     distances are measured in pixels. 
#     OUTPUT: 
#         tabout:
#             tabout: a table with columns: isocut, iso_d, iso_r, iso_area, iso_flux

#         iso_d:
#              the largest distance between any pair of superpixel that has flux > isocut.
#         iso_r:
#              the largest distance between image center and superpixels that has flux > isocut.
#         iso_area:
#              the sum of area in the superpixels that has flux > isocut.
#         iso_flux:
#              the sum of flux in the superpixels that has flux > isocut. 
#              (This may not capture all the faint emission, but will not 
#              suffer from negative over subtractions)

#     """
#     #==== protect the original params from being altered
#     from copy import deepcopy
#     img=deepcopy(img)
#     isocut=deepcopy(isocut)

#     #==== preprocessing
#     # sanity check for unit consistency
#     if useunits:
#         for var in [img,isocut]:
#             if hasattr(var, 'unit'):
#                 if var.unit != imgunit: raise ValueError('unit inconsistent')
#     else: 
#         pixelsize,pixunit,imgunit=1.,1.,1.

#     # if img has a unit get rid of it, for the convenience of calculation
#     if hasattr(img, 'unit'): img=img.value
#     if hasattr(isocut, 'unit'): isocut=isocut.value

#     #==== clipping
#     # set all the faint featuers 0
#     img[img<isocut]=0.
#     nx,ny=img.shape
#     nbp=np.sum(img>0.) # count number of bright super pixels

#     #==== calculation based on pure pixel values
#     if nbp ==0:
#         iso_d, iso_r, iso_area, iso_flux=0., 0., 0., 0.        
#     elif nbp >0:
#         # calculate iso_d
#         ixs,iys=np.where(img>0)
#         distances=np.zeros(nbp*nbp)
#         for i in range(nbp):
#             for j in range(nbp):
#                 distances[i*nbp+j]=np.sqrt((ixs[i]-ixs[j])**2+(iys[i]-iys[j])**2)
#         iso_d=max(distances)
#         # calculate iso_r
#         distances=np.zeros(nbp)
#         for i in range(nbp):
#             distances[i]=np.sqrt((ixs[i]-0.5*(nx-1))**2+(iys[i]-0.5*(ny-1))**2)
#         iso_r=max(distances)
#         # calculate iso_area
#         iso_area=nbp
#         # calculate iso_flux
#         iso_flux=np.sum(img)

#     # adding back the units
#     iso_d=iso_d*pixelsize*pixunit
#     iso_r=iso_r*pixelsize*pixunit
#     iso_area=iso_area*(pixelsize*pixunit)**2
#     iso_flux=iso_flux*imgunit*(pixelsize*pixunit)**2
#     isocut=isocut*imgunit

#     # output
#     tabout=Table()
#     for var in ['isocut','iso_d','iso_r','iso_area','iso_flux']: 
#         if isinstance(eval(var),u.quantity.Quantity):
#             tabout[var]=[eval(var).value]
#             tabout[var].unit=eval(var).unit
#         else: 
#             tabout[var]=[eval(var)]

#     return tabout
