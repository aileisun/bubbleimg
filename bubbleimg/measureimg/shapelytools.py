# shapelytools.py
"""
use shapely for contour differencing
"""
import shapely.geometry as sg
import numpy as np

import polytools

def diffcontoursdict(contoursdict1, contoursdict2, filecontoursdict=None):
    """
    return the contoursdict of contoursdict1 minus contoursdict2
    """
    import copy

    contours1=contoursdict1['contours']

    if polytools.NetPolygonsArea(contours1)>0.: 
        contours2=contoursdict2['contours']        
        contoursdiff=diffcontours(contours1, contours2)
    else: 
        contoursdiff=[]

    contoursdictout=copy.copy(contoursdict1)
    if filecontoursdict is not None:
        contoursdictout['filecontoursdict']=filecontoursdict
    contoursdictout['contours']=contoursdiff
    contoursdictout['masked']=True
    contoursdictout['mask']=contoursdict2['filecontoursdict']
    return contoursdictout

def diffcontours(contours1, contours2):
    """
    return the contours of contours1 minus contours2
    """
    sgM1=sgMultipolygon_from_contours(contours1)
    sgM2=sgMultipolygon_from_contours(contours2)

    sgMdiff=sgM1.difference(sgM2)

    if isinstance(sgMdiff, sg.polygon.Polygon):
        contours= contours_from_sgPolygon(sgMdiff)
    elif isinstance(sgMdiff, sg.MultiPolygon):
        contours=contours_from_sgMultipolygon(sgMdiff)
    elif isinstance(sgMdiff, sg.collection.GeometryCollection): 
        contours=contours_from_sgGeoCollection(sgMdiff)
    elif sgMdiff.area <=0:
        contours=[]
    else: 
        print type(sgMdiff)
        raise NameError('sgMdiff type not understood' )
    return contours

def contours_from_sgPolygon(sgP):
    """
    transform a sg polygon instance into contours
    """
    if not isinstance(sgP, sg.polygon.Polygon): 
        raise ValueError("input is not shapely Polygon")

    highcontour=np.array(sgP.exterior.coords)
    lowcontours=[np.array(interior.coords) for interior in sgP.interiors]

    # make sure the rotation sense is correct
    if polytools.SignedPolygonArea(highcontour)<0.: 
        highcontour=highcontour[::-1,:]
    for i in range(len(lowcontours)):
        if polytools.SignedPolygonArea(lowcontours[i])>0.: 
            lowcontours[i]=lowcontours[i][::-1,:]

    contours=[highcontour]+lowcontours
    return contours

def contours_from_sgMultipolygon(sgM):
    """
    transform a sg multi polygon instance into contours
    """
    if not isinstance(sgM, sg.MultiPolygon): 
        raise ValueError("input is not shapely MultiPolygon")
    contours=[]
    for sgP in sgM:
        contours = contours + contours_from_sgPolygon(sgP)
    return contours


def contours_from_sgGeoCollection(sgC):
    """ transform a sg geocollection instance into contours """
    if not isinstance(sgC, sg.collection.GeometryCollection): 
        raise ValueError("input is not shapely GeometryCollection")
    contours=[]
    for sgP in sgC:
        if isinstance(sgP, sg.polygon.Polygon):
            contours = contours + contours_from_sgPolygon(sgP)
    return contours



def sgMultipolygon_from_contours(contours):
    """
    transform contours to sg.MultiPolygon object

    assuming all the low contours are enclosed by high contours, 
    and high contours are not intersecting. 
    """
    highcontours=polytools.select_highcontours(contours)
    lowcontours=polytools.select_lowcontours(contours)

    ps=[]
    for highcontour in highcontours:
        holes=polytools.select_enclosedcontours(lowcontours, [highcontour])
        sgP=sg.Polygon(highcontour,holes)
        ps=ps+[sgP]

    sgM = sg.MultiPolygon(ps)
    return sgM



