# subtractexp.py
# ALS 2017/06/12

import copy

def subtract_exp_with_ratio(exp1, exp2, a1=1., a2=1.):
	"""
	calculate a1 * exp1 - a2 * exp2, where exp are lsst ExposureF. 
	The image, mask, and variance planes contains the new subtracted image. 
	The image information (higher extensions) is the same as exp1. 

	Params
	------
	exp1 (ExposureF)
	exp2 (ExposureF)
	a1=1. (float)
	a2=1. (float)

	Return
	------
	expout (ExposureF)
	"""
	exp1 = copy.copy(exp1)
	exp2 = copy.copy(exp2)

	mImg1 = exp1.getMaskedImage()
	mImg2 = exp2.getMaskedImage()

	mImg1 *= a1
	mImg2 *= a2

	mImg1 -= mImg2

	return exp1