# test_subtractexp.py


# import pytest
# import shutil
# import os
import numpy as np

import lsst.afw.image as afwImage
# from lsst.afw.image import ExposureF
from lsst.afw.geom import degrees
from lsst.afw.coord import IcrsCoord

from .. import subtractexp


ra = 140.099341430207
dec = 0.580162492432517
# dir_parent = './testing/'
# dir_obj = './testing/SDSSJ0920+0034/'

dir_verif = 'verification_data/SDSSJ0920+0034/'
# bandline = 'i'
# bandconti = 'r'
# survey = 'hsc'
# z = 0.4114188

def test_subtractexp():

	a1 = 0.5
	a2 = 2.

	exp1 = afwImage.ExposureF(dir_verif+'stamp-i.fits')
	exp2 = afwImage.ExposureF(dir_verif+'stamp-y.fits')

	img1 = exp1.getMaskedImage().getImage().getArray() # np array
	img2 = exp2.getMaskedImage().getImage().getArray()

	expsub = subtractexp.subtract_exp_with_ratio(exp1, exp2, a1=a1, a2=a2)

	imgsub = expsub.getMaskedImage().getImage().getArray() # np array

	# assert exp1 is not altered
	assert np.all(img1 == exp1.getMaskedImage().getImage().getArray())

	# assert subtracted img is correct
	assert np.all(imgsub == (img1*a1 - img2*a2))

	# assert subtracted var correct
	var1 = exp1.getMaskedImage().getVariance().getArray() # np array
	var2 = exp2.getMaskedImage().getVariance().getArray()
	varsub = expsub.getMaskedImage().getVariance().getArray()

	assert np.all(varsub == (a1**2) * var1 + (a2**2) * var2)

	# assert filter of subtracted is the same as exp1
	assert expsub.getFilter().getName() == exp1.getFilter().getName()

	# assert psf corresponds to exp1
	assert np.all(get_psfImageKer(expsub) == get_psfImageKer(exp1))



def get_psfImageKer(exp):
	psf = exp.getPsf()
	w = exp.getWcs()
	c = IcrsCoord(ra*degrees, dec*degrees)
	pos = w.skyToPixel(c)

	psfImageKer = psf.computeKernelImage(pos)
