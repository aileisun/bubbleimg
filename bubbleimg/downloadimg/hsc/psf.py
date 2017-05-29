# psf.py
# ALS 2017/05/09

"""
extract psf from hsc image at location ra, dec. 
"""
from lsst.afw.image import ExposureF
from lsst.afw.geom import degrees
from lsst.afw.coord import IcrsCoord


def exposureF_to_psf(fn_in, fn_out, ra, dec ):
	"""
	extract psf from hsc image at location ra, dec, and write it to file "fn_out".

	Params
	------
	fn_in (str): filepath towards an ExposureF file. 
	fn_out (str)
	ra (float)
	dec (float)

	"""

	exp = ExposureF(fn_in)
	psf = exp.getPsf()

	w = exp.getWcs()
	c = IcrsCoord(ra*degrees, dec*degrees)
	pos = w.skyToPixel(c)

	psfImageKer = psf.computeKernelImage(pos)
	psfImageKer.writeFits(fn_out)

