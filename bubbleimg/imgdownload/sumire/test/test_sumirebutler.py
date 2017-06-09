##

import matplotlib
matplotlib.use('Agg')

import lsst.afw.image as afwImage
from lsst.afw.geom import degrees
from lsst.afw.coord import IcrsCoord
import os

# from lsst.daf.persistence import Butler
from ..sumirebutler import Butler

rerunpath = '/array2/SSP/dr1/s16a/data/s16a_wide/'

ra = 140.099364238908123
dec = 0.580160150759375104
tract = 9564
patch_s = "7,3"
# fn_in = "/array2/SSP/dr1/s16a/data/s16a_wide/deepCoadd/HSC-R/9564/7,3/calexp-HSC-R-9564-7,3.fits"

dataId = dict(tract=tract, patch=patch_s, filter="HSC-R")
# dataId = dict(tract=9564, patch="7,3", filter="HSC-R")


def test_sumirebutler_get():
	b = Butler(rerunpath)

	assert b.datasetExists("deepCoadd_calexp", dataId)

	exp = b.get("deepCoadd_calexp", dataId, immediate=True)
	assert isinstance(exp, afwImage.imageLib.ExposureF)


def test_sumirebutler_get_gzfiles():
	""" get file even if its compressed """
	assert False

	b = Butler(rerunpath)

	assert b.datasetExists("deepCoadd_calexp", dataId)

	exp = b.get("deepCoadd_calexp", dataId, immediate=True)
	assert isinstance(exp, afwImage.imageLib.ExposureF)
