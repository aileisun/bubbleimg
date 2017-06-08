# psf.py
# ALS 2017/05/09

"""
extract psf from hsc image at location ra, dec. 
"""
import os 
import gzip
import shutil
import sys

from lsst.afw.image import ExposureF
from lsst.afw.geom import degrees
from lsst.afw.coord import IcrsCoord

fn_temp = '/data/home/hscpipe/alsun/get_psf/temp/temp.fits'


def main(fn_in, fn_out, ra, dec ):
	"""
	extract psf from hsc image at location ra, dec, and write it to file "fn_out".

	Params
	------
	fn_in (str): filepath towards an ExposureF file. 
	fn_out (str)
	ra (float)
	dec (float)

	Return 
	------
	status (bool)

	"""

	ra = float(ra)
	dec = float(dec)

	if os.path.isfile(fn_out): 
		os.remove(fn_out)

	# deal with case when the file is compressed as .gz -> uncompress and save as fn_temp
	if not os.path.isfile(fn_in):
		fn_in_gz = fn_in+'.gz'
		if os.path.isfile(fn_in_gz):
			if os.path.isfile(fn_temp):
				os.remove(fn_temp)
			ungz(fn_in_gz, fn_temp)
			fn_in_update = fn_temp
		else:
			print("False")
			return False
	else:
		fn_in_update = fn_in

	exp = ExposureF(fn_in_update)
	psf = exp.getPsf()

	w = exp.getWcs()
	c = IcrsCoord(ra*degrees, dec*degrees)
	pos = w.skyToPixel(c)

	psfImageKer = psf.computeKernelImage(pos)
	psfImageKer.writeFits(fn_out)

	status = os.path.isfile(fn_out)
	print(str(status))
	return status


def ungz(fn_in, fn_out):
	""" decompress .gz file """
	with gzip.open(fn_in, 'rb') as f_in, open(fn_out, 'wb') as f_out:
		shutil.copyfileobj(f_in, f_out)


if __name__ == '__main__':
	main(*sys.argv[1:])