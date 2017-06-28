# cutout.py
# ALS 2017/06/09
# from example code by Yusra AlSayyad

import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom

def make_cutout(exp, ra, dec, npix=128):
	""" 
	make cutout from ExposureF object and return the result as an ExposureF object 

	Params
	------
	exp (afwImage.ExposureF object)
	ra (float)
	dec (float)
	npix=128 (int)
	"""
	lsstwcs = exp.getWcs()
	pointCoord = afwCoord.IcrsCoord(afwGeom.Angle(ra, afwGeom.degrees), afwGeom.Angle(dec, afwGeom.degrees))
	x, y = lsstwcs.skyToPixel(pointCoord)
	corner = afwGeom.Point2I(int(x - npix // 2), int(y - npix // 2))
	bbox = afwGeom.Box2I(afwGeom.Point2I(corner), afwGeom.Extent2I(npix, npix))
	stamp = exp.Factory(exp, bbox, afwImage.PARENT, False)

	return stamp
