# matchpsf.py
# ALS 2017/06/10

def match_two_exp_psf(exp1, exp2):
	"""
	Return a pair of exposures that corresponds to exp1 and exp2 but have the same psf.
	Determine which one has smaller psf and match that to the other. 

	Params
	------
	exp1 (ExposureF)
	exp2 (ExposureF)

	Return
	------
	resultExp1 (ExposureF)
	resultExp2 (ExposureF)
	"""
	if has_smaller_psf(exp1, exp2):
		expm1 = match_psf(exp1, exp2)
		return expm1, exp2
	elif has_smaller_psf(exp2, exp1):
		expm2 = match_psf(exp2, exp1)
		return exp1, expm2
	else:
		print("[matchpsf] skip matching psf as the psfs are similar")
		return exp1, exp2


def match_psf(fromExp, toExp): 
	"""
	match the psf of the fromExp to that of toExp
	return resultExp that has image convolved with kernel and psf corresponding to the convolved psf. 

	Params
	------
	fromExp
	toExp

	Return
	------
	resultExp
	"""

	raise NotImplementedError("to be written")

 #    config = ImagePsfMatchTask.ConfigClass()
 #    config.kernel.name = "AL"
 #    config.kernel.active.fitForBackground = True
 #    config.kernel.active.spatialKernelOrder = 1
 #    config.kernel.active.spatialBgOrder = 0
 #    config.kernel.active.sizeCellX = 128
 #    config.kernel.active.sizeCellY = 128

	# psfMatchTask = ImagePsfMatchTask(config=config)
	# result = psfMatchTask.matchExposures(templateExp, scienceExp)

	# result = psfMatchTask.subtractExposures(templateExp, scienceExp)

	return exp

def has_smaller_psf(exp1, exp2, diffInPix=0.01):
	""" whether exp1 has psf smaller than exp2 by a sigma more than diffInPix """

	sig1 = getPsfSigmaInPix(exp1.getPsf())
	sig2 = getPsfSigmaInPix(exp2.getPsf())

	result = sig1 < sig2 - diffInPix

	return result



def getPsfSigmaInPix(psf):
    """ Return the sigma in pixels of a Psf """
    sigPix = psf.computeShape().getDeterminantRadius()
    return sigPix

