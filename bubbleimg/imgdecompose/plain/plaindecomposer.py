# decomposer.py 
# ALS 2017/06/01
import os
from astropy.io import fits

from .. import Decomposer
import matchpsf

class plainDecomposer(Decomposer):

	def __init__(self, **kwargs):
		"""
		child of decomposer

		Use homemade code to do psf matching. Requires psf-{band}.fits files. 
		"""
		
		super(plainDecomposer, self).__init__(**kwargs)


	def get_fp_psf_kernelcnvl(self, band, bandto):
		return self.dir_obj+'psf-{0}_kernelcnvl-{1}.fits'.format(band, bandto)


	def get_fp_psf_psfmatched(self, band, bandto):
		return self.dir_obj+'psf-{0}_psfmatched-{1}.fits'.format(band, bandto)


	def make_stamp_contsub(self, band, bandconti, overwrite=True):
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
			print("[decomposer] making stamp_contsub")

			ratioconti = self._get_conti_fnu_ratio_from_spector(band, bandconti)

			if self._has_smaller_psf(band, bandconti):
				print("[decomposer] matching psf of {}-band to conti {}-band".format(band, bandconti))
				self._make_stamp_psfmatch(band=band, bandto=bandconti, overwrite=overwrite)
				fn_band = self.get_fp_stamp_psfmatched(band=band, bandto=bandconti)
				fn_cont = self.get_fp_stamp(bandconti)
			elif self._has_smaller_psf(bandconti, band):
				print("[decomposer] matching psf of conti {}-band to {}-band".format(bandconti, band))
				self._make_stamp_psfmatch(band=bandconti, bandto=band, overwrite=overwrite)
				fn_band = self.get_fp_stamp(band)
				fn_cont = self.get_fp_stamp_psfmatched(band=bandconti, bandto=band)
			else: 
				print("[decomposer] skip matching psf as they are similar")
				fn_band = self.get_fp_stamp(band)
				fn_cont = self.get_fp_stamp(bandconti)

			self._subtract_img_w_ratio(fn1=fn_band, fn2=fn_cont, fnout=fn, a1=1., a2=ratioconti, overwrite=overwrite)

		else:
			print("[decomposer] skip making stamp_contsub as files exist")

		return os.path.isfile(fn)


	def _make_stamp_psfmatch(self, band, bandto, overwrite=True):
		""" 
		make file stamp-x_psfmatch-y.fits -- img stamp-x matched to the psf of y.
		Find best fit moffat convolving kernel which when convolved with psf-x becomes psf-y. 

		Params
		------
		self
		band (str)
		bandto (str)
		overwrite=False

		Return
		------
		status

		Write Output (e.g., if band = 'i', bandto = 'z')
		------------
		stamp-i_psfmatched-z.fits
		psf-i_cnvlkernel-z.fits
		psf-i_psfmatched-z.fits
		"""

		fns = [self.get_fp_stamp_psfmatched(band, bandto), 
				self.get_fp_psf_psfmatched(band, bandto),
				self.get_fp_psf_kernelcnvl(band, bandto),
				]

		if not all([os.path.isfile(fn) for fn in fns]) or overwrite:
			print("[decomposer] making stamp_psfmatch")
			psf = fits.getdata(self.get_fp_psf(band))
			psfto = fits.getdata(self.get_fp_psf(bandto))

			hdus = fits.open(self.get_fp_stamp(band))
			stamp = hdus[0].data

			stamp_cnvled, psf_cnvled, kernel_cnvl = matchpsf.match_psf(stamp, psf, psfto)

			# write psf matched stamp
			hdus[0].data = stamp_cnvled
			hdus[0].header['BDPSFTO'] = (bandto, 'the band the PSF is matched to')
			hdus[0].header['COMMENT'] = "PSF matched by ALS"
			hdus.writeto(self.get_fp_stamp_psfmatched(band, bandto), overwrite=overwrite)

			# write psf matched psf
			fits.PrimaryHDU(psf_cnvled).writeto(self.get_fp_psf_psfmatched(band, bandto), overwrite=overwrite)

			# write convolving kernel
			fits.PrimaryHDU(kernel_cnvl).writeto(self.get_fp_psf_kernelcnvl(band, bandto), overwrite=overwrite)
		else:
			print("[decomposer] skip making stamp_psfmatch as files exist")

		return all([os.path.isfile(fn) for fn in fns])



	def _has_smaller_psf(self, band, bandto, diffpsf_threshold=0.01):
		"""
		whether band has smaller psf than bandto by a margin of "diffpsf_threshold" in arcsec. 
		get psf size from hsc_xid if the survey = 'hsc'. 

		Params
		------
		band (str)
		bandto (str)
		diffpsf_threshold=0.1:
			the threshold of psf difference in arcsec. If the diff is larger than return True. 

		Return
		------
		answer (bool)
		"""

		if self.survey == 'hsc':
			psize_from = self.obj.hsc.get_psfsize(band=band)
			psize_to = self.obj.hsc.get_psfsize(band=bandto)

			return psize_to - psize_from > diffpsf_threshold
		else:
			raise Exception("[decomposer] _has_smaller_psf() for this survey is not constructed")


	def _subtract_img_w_ratio(self, fn1, fn2, fnout, a1=1., a2=1., overwrite=False):
 
		if not os.path.isfile(fnout) or overwrite:
			hdu1 = fits.open(fn1)
			img1 = hdu1[0].data

			img2 = fits.getdata(fn2)

			imgout = a1*img1-a2*img2

			hduout = hdu1
			hduout[0].data = imgout
			hduout[0].header['COMMENT'] = 'image subtracted by ALS'
			hduout[0].header['COMMENT'] = '{}*{} minus {}*{}'.format('%.2f'%a1, fn1.split('/')[-1], '%.2f'%a2, fn2.split('/')[-1])
			hduout.writeto(fnout, overwrite=overwrite)
