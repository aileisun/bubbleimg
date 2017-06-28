# decomposer.py 
# ALS 2017/06/01
import os
from astropy.io import fits
import lsst.afw.image as afwImage

from .. import Decomposer

import matchpsf
import subtractexp

class lsstDecomposer(Decomposer):

	def __init__(self, **kwargs):
		"""
		child of decomposer

		Use homemade code to do psf matching. Requires psf-{band}.fits files. 
		"""
		
		super(lsstDecomposer, self).__init__(**kwargs)


	def _get_stamp_exposure(self, band):
		fp = self.get_fp_stamp(band)
		return afwImage.ExposureF(fp)


	def _get_contsub_exposure(self, band, bandconti):
		fp = self.get_fp_stamp_contsub(band, bandconti)
		return afwImage.ExposureF(fp)


	def make_stamp_linemap(self, bandline, bandconti, line='OIII5008', overwrite=False):
		"""
		make stamp of line map in observed frame (?) flux with units [erg s-1 cm-2]

		Params
		------
		self
		bandline  (str)
		bandconti (str)
		line = 'OIII5008' (str)
		overwrite = False (bool)

		Return
		------
		status (bool)

		Write Output 
		------------
		e.g., stamp-lOIII5008.fits
		"""
		fn = self.get_fp_stamp_line(line)

		if not os.path.isfile(fn) or overwrite:
			self.make_stamp_contsub(band=bandline, bandconti=bandconti, overwrite=overwrite)

			s = self._get_spector()
			r = s.calc_fline_over_fnuband(band=bandline, line=line)

			exp = self._get_contsub_exposure(band=bandline, bandconti=bandconti)
			img = exp.getMaskedImage()
			img *= r

			raise Exception("need to check whether ratio with unit can be multipled by maskedImg")

			exp.writeFits(fn)

		else:
			print("[lsstdecomposer] skip making stamp_contsub as files exist")

		return os.path.isfile(fn)


	def make_stamp_contsub(self, band, bandconti, overwrite=False):
		"""
		make stamp that is continuum subtracted

		Params
		------
		self
		band (str)
		bandconti (str)
		overwrite=False

		Return
		------
		status

		Write Output (e.g., if band = 'i', bandto = 'z')
		------------
		stamp-i_contsub-z.fits
		"""
		fn = self.get_fp_stamp_contsub(band, bandconti)

		if not os.path.isfile(fn) or overwrite:
			print("[lsstdecomposer] making stamp_contsub")

			ratioconti = self._get_conti_fnu_ratio_from_spector(band, bandconti)

			exp = self._get_stamp_exposure(band)
			exp_conti = self._get_stamp_exposure(bandconti)

			expm, expm_conti = matchpsf.match_two_exp_psf(exp, exp_conti)
			exp_sub = subtractexp.subtract_exp_with_ratio(expm, expm_conti, a1=1., a2=ratioconti)

			exp_sub.writeFits(fn) 

		else:
			print("[lsstdecomposer] skip making stamp_contsub as files exist")

		return os.path.isfile(fn)



		
	# def make_stamp_psfmatch(self, band, bandto, overwrite=True):
	# 	""" 
	# 	make stamp that has psf matched to stamp of another band

	# 	Params
	# 	------
	# 	self
	# 	band (str)
	# 	bandto (str)
	# 	overwrite=False

	# 	Return
	# 	------
	# 	status

	# 	Write Output (e.g., if band = 'i', bandto = 'z')
	# 	------------
	# 	stamp-i_psfmatched-z.fits
	# 	"""
	# 	fn = self.get_fp_stamp_psfmatched(band, bandto)

	# 	if not os.path.isfile(fn) or overwrite:
	# 		print("[lsstdecomposer] making stamp_psfmatch")
	# 		raise NotImplementedError("to be written")

	# 		# psf = fits.getdata(self.get_fp_psf(band))
	# 		# psfto = fits.getdata(self.get_fp_psf(bandto))

	# 		# hdus = fits.open(self.get_fp_stamp(band))
	# 		# stamp = hdus[0].data

	# 		# stamp_cnvled, psf_cnvled, kernel_cnvl = matchpsf.match_psf(stamp, psf, psfto)

	# 		# # write psf matched stamp
	# 		# hdus[0].data = stamp_cnvled
	# 		# hdus[0].header['BDPSFTO'] = (bandto, 'the band the PSF is matched to')
	# 		# hdus[0].header['COMMENT'] = "PSF matched by ALS"
	# 		# hdus.writeto(self.get_fp_stamp_psfmatched(band, bandto), overwrite=overwrite)

	# 		# # write psf matched psf
	# 		# fits.PrimaryHDU(psf_cnvled).writeto(self.get_fp_psf_psfmatched(band, bandto), overwrite=overwrite)

	# 		# # write convolving kernel
	# 		# fits.PrimaryHDU(kernel_cnvl).writeto(self.get_fp_psf_kernelcnvl(band, bandto), overwrite=overwrite)
	# 	else:
	# 		print("[lsstdecomposer] skip making stamp_psfmatch as files exist")

	# 	return os.path.isfile(fn)




