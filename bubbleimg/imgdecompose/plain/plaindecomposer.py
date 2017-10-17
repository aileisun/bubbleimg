# decomposer.py 
# ALS 2017/06/01
import os
from shutil import copyfile
from astropy.io import fits
import astropy.units as u
import astropy.table as at
import numpy as np

from ..decomposer import Decomposer
from ... import tabtools

import matchpsf
import plotpsf_tools

class plainDecomposer(Decomposer):

	def __init__(self, **kwargs):
		"""
		child of decomposer

		Use homemade code to do psf matching. Requires psf-{band}.fits files. 
		"""
		
		super(plainDecomposer, self).__init__(**kwargs)


	def get_fp_psf_tab(self, msrsuffix=''):
		return self.dir_obj + 'psf{}.csv'.format(msrsuffix)


	def get_fp_psf(self, fp_stamp):
		""" return the path of a psf file corresponding to that of the stamp """

		dirpath = os.path.dirname(fp_stamp)
		fn_stamp = os.path.basename(fp_stamp)

		fn_psf = 'psf'+fn_stamp[5:]
		fp_psf = os.path.join(dirpath, fn_psf)

		return fp_psf


	def get_fp_psk(self, fp_stamp):
		""" return the path of a psf matching kernel file corresponding to that of the stamp """

		dirpath = os.path.dirname(fp_stamp)
		fn_stamp = os.path.basename(fp_stamp)

		fn_psk = 'psk'+fn_stamp[5:]
		fp_psk = os.path.join(dirpath, fn_psk)

		return fp_psk


	def get_fp_contsubtab(self, band, bandconti):
		return self.dir_obj + 'contsub_{}{}.csv'.format(band, bandconti)


	def get_fp_contsubplot(self, band, bandconti):
		return self.dir_obj + 'contsub_psf_{}{}.pdf'.format(band, bandconti)


	def make_stamp_linemap_I(self, bandline, bandconti, line='OIII5008', overwrite=False):
		""" 
		make stamp of line map in rest frame intensity in units of [erg s-1 cm-2 arcsec-2]
		ready for isophotal measurements. 

		See make_stamp_linemap for details. 
		"""

		# setting
		fn_in = self.get_fp_stamp_line(line)
		fn_out = self.get_fp_stamp_line_I(line)

		if not os.path.isfile(fn_out) or overwrite:
			print("[decomposer] making stamp_linemap_I")

			self.make_stamp_linemap(bandline, bandconti, line=line, overwrite=True)
			self._convert_linemap_to_linemapI(fn_in, fn_out)

			self._copy_psf(fn_in, fn_out)
		else:
			print("[decomposer] skip making stamp_linemap_I as files exist")

		status = os.path.isfile(fn_out)
		return status


	def make_stamp_linemap(self, bandline, bandconti, line='OIII5008', overwrite=False):
		"""
		make stamp of line map in observed frame with flux in units of [erg s-1 cm-2]

		Params
		------
		self
		bandline  (str)
		bandconti (str)
		line = 'OIII5008' (str)
		overwrite = False (bool)
		mapUnit = u.Unit('1e-17 erg s-1 cm-2') (unit)

		Return
		------
		status (bool)

		Write Output 
		------------
		e.g., stamp-OIII5008.fits
		"""
		fn_in = self.get_fp_stamp_contsub(band=bandline, bandconti=bandconti)
		fn_out = self.get_fp_stamp_line(line)

		if not os.path.isfile(fn_out) or overwrite:
			print("[decomposer] making stamp_linemap")

			# prep file
			self.make_stamp_contsub(band=bandline, bandconti=bandconti, overwrite=False)
			self._convert_contsub_to_linemap(fn_in, fn_out, bandline, line)

			self._copy_psf(fn_in, fn_out)

		else:
			print("[decomposer] skip making stamp_linemap as files exist")

		status = os.path.isfile(fn_out)
		return status


	def make_stamp_contsub(self, band, bandconti, fwhm_fracdiff_threshold=0.00, toplot_psf=False, towrite_tab=False, overwrite=False):
		"""
		make stamp that is continuum subtracted

		Params
		------
		self
		band (str)
		bandconti (str)
		fwhm_fracdiff_threshold = 0.00
			the psf fwhm fractional difference above which psf matching will be performed
		toplot_psf=False
			whether to make plots psf-*_psfmt-*.pdf showing the psf matching performance
		towrite_tab = False
			whether to make tab contsub-{}-{}.csv which summarizes which band is psf matched to which band
		overwrite=False

		Return
		------
		status

		Write Output (e.g., if band = 'i', bandto = 'z')
		------------
		stamp-i_contsub-z.fits
		"""
		fn_out = self.get_fp_stamp_contsub(band, bandconti)

		if not os.path.isfile(fn_out) or overwrite:
			print("[decomposer] making stamp_contsub")

			ratioconti = self._get_conti_fnu_ratio_from_spector(band, bandconti)

			# psf matching accorindg to which psf is larger
			if self._band_has_smaller_psf(band, bandconti, fracdiff_threshold=fwhm_fracdiff_threshold):
				print("[decomposer] matching psf of {}-band to conti {}-band".format(band, bandconti))
				bandfrom = band
				bandto = bandconti
				fn_band = self.get_fp_stamp_psfmatched(band=bandfrom, bandto=bandto)
				fn_cont = self.get_fp_stamp(bandto)

			elif self._band_has_smaller_psf(bandconti, band, fracdiff_threshold=fwhm_fracdiff_threshold):
				print("[decomposer] matching psf of conti {}-band to {}-band".format(bandconti, band))
				bandfrom = bandconti
				bandto = band
				fn_band = self.get_fp_stamp(bandto)
				fn_cont = self.get_fp_stamp_psfmatched(band=bandfrom, bandto=bandto)

			else: 
				print("[decomposer] skip matching psf as they are similar")
				bandfrom = None
				bandto = None
				fn_band = self.get_fp_stamp(band)
				fn_cont = self.get_fp_stamp(bandconti)

			# operation
			self.make_stamp_psfmatch(band=bandfrom, bandto=bandto, overwrite=overwrite)
			self._subtract_img_w_ratio(fn1=fn_band, fn2=fn_cont, fnout=fn_out, a1=1., a2=ratioconti, overwrite=overwrite)



			# copy psf
			fn_psf_in = self.get_fp_psf(fn_band)
			fn_psf_out = self.get_fp_psf(fn_out)
			copyfile(fn_psf_in, fn_psf_out)

			if toplot_psf:
				fn_plot = self.get_fp_contsubplot(band, bandconti)
				if (bandfrom is None) or (bandto is None):
					self.plot_psfmatch(band=band, bandto=bandconti, fn=fn_plot, matching=False, overwrite=overwrite)
				else: 
					self.plot_psfmatch(band=bandfrom, bandto=bandto, fn=fn_plot, matching=True, overwrite=overwrite)

			if towrite_tab:
				fn_tab = self.get_fp_contsubtab(band, bandconti)
				tab = self._get_contsub_tab_content(band, bandconti, bandfrom, bandto, fwhm_fracdiff_threshold)
				tab.write(fn_tab, overwrite=True)

		else:
			print("[decomposer] skip making stamp_contsub as files exist")

		return os.path.isfile(fn_out)


	def make_psf_tab_bands(self, mode='moffat', overwrite=False):
		"""
		make measurements on the psf sizes of stamp-*.fits of all bands and write to file 'psf.csv'. 
		"""
		fn = self.get_fp_psf_tab(msrsuffix='_bands')

		if not os.path.isfile(fn) or overwrite:
			print("[plaindecomposer] making psf.csv")
			tab = at.Table([[mode]], names=['mode'])

			for band in self.bands: 
				psf_fwhm_arcs = self._get_psf_fwhm_arcs(imgtag=band, mode=mode)
				tab['psf_fwhm_arcs_{}'.format(band)] = psf_fwhm_arcs

			for band in self.bands: 
				psf_fwhm_pix = self._get_psf_fwhm_pix(imgtag=band, mode=mode)
				tab['psf_fwhm_pix_{}'.format(band)] = psf_fwhm_pix

			tab.write(fn, format='ascii.csv', overwrite=overwrite)
		else:
			print("[plaindecomposer] skip making psf.csv as file exists")

		status = os.path.isfile(fn)
		return status


	def make_psf_tab(self, imgtag='OIII5008_I', mode='moffat', msrsuffix='', overwrite=False, append=False):
		"""
		make measurements on the psf sizes and write to file "psf{msrsuffix}.csv". 

		Params
		------
		imgtag='OIII5008_I'
			indicate which image to measure on
		mode='moffat'
		msrsuffix=''
			the suffix for the measurement file name
		overwrite=False
		append=False
			if true than always append line to end of file
		"""
		fn = self.get_fp_psf_tab(msrsuffix=msrsuffix)
		condi = dict(imgtag=imgtag, mode=mode)

		if append or overwrite or (not tabtools.fn_has_row(fn, condi)):
			tab = at.Table([[imgtag], [mode]], names=['imgtag', 'mode'])

			psf_fwhm_arcs = self._get_psf_fwhm_arcs(imgtag=imgtag, mode=mode)
			tab['psf_fwhm_arcs'] = psf_fwhm_arcs

			psf_fwhm_pix = self._get_psf_fwhm_pix(imgtag=imgtag, mode=mode)
			tab['psf_fwhm_pix'] = psf_fwhm_pix

			if mode == 'moffat':
				gamma, alpha = self._get_psf_moffat_params(imgtag)
				tab['moffat_gamma'] = gamma
				tab['moffat_alpha'] = alpha

			tabtools.write_row(fn=fn, row=tab, condi=condi, overwrite=overwrite, append=append)

		else:
			print("[plaindecomposer] skip making psf measurement as files exist")

		return os.path.isfile(fn)


	def make_stamp_psfmatch(self, band, bandto, overwrite=False, towrite_psk=False):
		""" 
		make file stamp-x_psfmatch-y.fits -- img stamp-x matched to the psf of y.
		Find best fit moffat convolving kernel which when convolved with psf-x becomes psf-y. 

		Params
		------
		self
		band (str)
		bandto (str)
		overwrite=False (bool)
		towrite_psk=False (bool)
			whether to write psf matching kernel to file

		Return
		------
		status (bool)

		Write Output (e.g., if band = 'i', bandto = 'z')
		------------
		stamp-i_psfmt-z.fits
		psf-i_psfmt-z.fits
		psk-i_psfmt-z.fits
		"""

		if (band is not None) and (bandto is not None):
			# define file names
			fp_img = self.get_fp_stamp(band)  # input stamp
			fp_psf = self.get_fp_psf(fp_img)  # input psf
			fp_psfto = self.get_fp_psf(self.get_fp_stamp(bandto))  # psf to match to

			fp_img_out = self.get_fp_stamp_psfmatched(band, bandto)
			fp_psf_out = self.get_fp_psf(fp_img_out)
			fp_psk_out = self.get_fp_psk(fp_img_out)

			# operation
			if not os.path.isfile(fp_img_out) or overwrite:
				matchpsf.match_psf_fits(fp_img, fp_psf, fp_psfto, fp_img_out, fp_psf_out, fp_psk_out, overwrite=overwrite, towrite_psk=towrite_psk)
			else:
				print("[decomposer] skip making stamp_psfmatch as files exist")

			status = all([os.path.isfile(fn) for fn in [fp_img_out, fp_psf_out]])
		else:
			print("[decomposer] skip make_stamp_psfmatch as input bands is None")
			status = True

		return status


	def plot_psfmatch(self, band, bandto, fn, matching, overwrite=False):
		""" make plot to visualize psf matching """
		if not os.path.isfile(fn) or overwrite:
			psfs = self._get_normalized_trimmed_psfs_for_plot(band=band, bandto=bandto, matching=matching)

			if psfs == None:
				print("[decomposer] plotting psf failed as files do not exist")
				status = False

			else: 
				print("[decomposer] plotting psf")
				(psf0, psf1, psfmt) = psfs

				plotpsf_tools._plot_psfmatch_kernel(band, bandto, psf0, psf1, psfmt, fn)
				status = os.path.isfile(fn)
		else:
			print("[decomposer] skip plotting psf as files exist")
			status = True

		return status


	def _get_contsub_tab_content(self, band, bandconti, bandfrom, bandto, fwhm_fracdiff_threshold, mode='moffat'):
		"""
		return the table that contains info about the contsub prosedure 
		"""
		psf_fwhm = self._get_psf_fwhm_arcs(imgtag=band, mode=mode)
		psf_fwhm_conti = self._get_psf_fwhm_arcs(imgtag=bandconti, mode=mode)

		fwhm_fracdiff = np.absolute(psf_fwhm - psf_fwhm_conti)/max(psf_fwhm, psf_fwhm_conti)

		if (bandfrom is not None) and (bandto is not None):

			psf_matching_err = self._calc_psf_matching_err(bandfrom=bandfrom, bandto=bandto, matching=True)
			tab = at.Table([[mode], [fwhm_fracdiff_threshold], [band], [bandconti], [psf_fwhm], [psf_fwhm_conti], [fwhm_fracdiff], [bandfrom], [bandto], [psf_matching_err]], names=['mode', 'fwhm_fracdiff_threshold', 'band', 'bandconti', 'psf_fwhm', 'psf_fwhm_conti', 'fwhm_fracdiff', 'band_psfm_from', 'band_psfm_to', 'psf_matching_err'])
		else:
			psf_matching_err = self._calc_psf_matching_err(bandfrom=band, bandto=bandconti, matching=False)
			tab = at.Table([[mode], [fwhm_fracdiff_threshold], [band], [bandconti], [psf_fwhm], [psf_fwhm_conti], [fwhm_fracdiff], ['None'], ['None'], [psf_matching_err]], names=['mode', 'fwhm_fracdiff_threshold', 'band', 'bandconti', 'psf_fwhm', 'psf_fwhm_conti', 'fwhm_fracdiff', 'band_psfm_from', 'band_psfm_to', 'psf_matching_err'])

		return tab


	def _calc_psf_matching_err(self, bandfrom, bandto, matching=True):
		"""
		calculate the maximum difference between the original psf and the one matched to it (psfmt) shown as a fraction of the maximum of the normalized psfs
		"""
		if matching:
			fp_psf_from = self.get_fp_psf(self.get_fp_stamp_psfmatched(bandfrom, bandto))
		else: 
			fp_psf_from = self.get_fp_psf(self.get_fp_stamp(bandfrom))

		fp_psf_to = self.get_fp_psf(self.get_fp_stamp(bandto))

		frac_diff = matchpsf.calc_max_frac_diff_between_two_psf(fp_psf_from, fp_psf_to)
		return frac_diff


	def _get_normalized_trimmed_psfs_for_plot(self, band, bandto, matching=True):
		try:
			psf0 = fits.getdata(self.get_fp_psf(self.get_fp_stamp(band=band)))
			psf1 = fits.getdata(self.get_fp_psf(self.get_fp_stamp(band=bandto)))
			if matching:
				psfmt = fits.getdata(self.get_fp_psf(self.get_fp_stamp_psfmatched(band, bandto)))
			else:
				psfmt = psf0

		except IOError:
			return None

		else: 
			psf0, psf1, psfmt = plotpsf_tools._normalize_trim_psf(psf0, psf1, psfmt)
			return (psf0, psf1, psfmt)


	def _convert_linemap_to_linemapI(self, fn_in, fn_out, mapUnit=u.Unit('1e-15 erg s-1 cm-2 arcsec-2')):
		"""
		take stamp_linemap in observed frame flux and convert it to stamp_linemap in rest frame intensity
		write files. Takes self.z. 

		uses F_obs/omega = I_obs    --- *(1+z)^-4  --->  I_restframe
		always overwrites

		Params
		------
		fn_in  (str): 
			input file path to linemap obs frame flux fits file
		fn_out (str):
			output file path to linemap rest frame intensity fits file
		mapUnit=u.Unit('1e-15 erg s-1 cm-2 arcsec-2') (unit)

		"""

		# read in
		hdus = fits.open(fn_in)
		img_in = hdus[0].data * u.Unit(hdus[0].header['BUNIT'])

		# calc ratio
		invOmega = 1./(self.pixsize**2)
		zdimming = (1.+self.z)**4

		# scale to line flux
		img_out = invOmega*zdimming*img_in
		img_out = img_out.to(mapUnit)

		# prepare new hdu
		hdus[0].data = img_out.value
		bunit = img_out.unit.to_string()
		hdus[0].header.set(keyword='BUNIT', value=bunit, comment="in rest frame")
		hdus[0].header.set(keyword='COMMENT', value="Scaled to rest frame line intensity by ALS")
		hdus[0].header.set(keyword='COMMENT', value="    assuming redshift {}".format('%.3f'%self.z))

		# write to fits file
		hdus.writeto(fn_out, overwrite=True)


	def _convert_contsub_to_linemap(self, fn_in, fn_out, bandline, line, mapUnit=u.Unit('1e-17 erg s-1 cm-2')):
		"""
		scale contsub to linemap (flux in obs frame) of the specific line taking into account the fraction of flux in different lines. Assuming all (important) lines have similar spatial distribution. See spector.calc_fline_over_fnuband() for more details. 

		always overwrites

		Params
		------
		fn_in  (str): 
			input file path to stamp_contsub
		fn_out (str):
			output file path to stamp_linemap 
		bandline (str)
		line (str)
		mapUnit=u.Unit('1e-17 erg s-1 cm-2') (unit)

		"""

		# read in file
		hdus = fits.open(fn_in)
		img_band = hdus[0].data * u.Unit(hdus[0].header['BUNIT'])

		# calc ratio
		s = self._get_spector()
		r_fline_over_fnuband = s.calc_fline_over_fnuband(band=bandline, line=line)

		# scale to line flux
		img_line = r_fline_over_fnuband	* img_band	
		img_line = img_line.to(mapUnit)

		# prepare new hdu
		hdus[0].data = img_line.value
		bunit = img_line.unit.to_string()
		hdus[0].header.set(keyword='BUNIT', value=bunit, comment="in observed frame")
		hdus[0].header.set(keyword='COMMENT', value="Scaled to {} line flux by ALS".format(line))
		hdus[0].header.set(keyword='COMMENT', value="    with ratio {} applied to {} band".format(str(r_fline_over_fnuband), bandline))
		hdus[0].header.set(keyword='COMMENT', value="    bunit converted from nanomaggy to {} in obs frame".format(bunit))

		# write to fits file
		hdus.writeto(fn_out, overwrite=True)


	def _subtract_img_w_ratio(self, fn1, fn2, fnout, a1=1., a2=1., overwrite=False):
 
		if not os.path.isfile(fnout) or overwrite:
			hdu1 = fits.open(fn1)
			img1 = hdu1[0].data

			img2 = fits.getdata(fn2)

			imgout = a1*img1-a2*img2

			hduout = hdu1
			hduout[0].data = imgout
			hduout[0].header['COMMENT'] = 'Image subtracted by ALS'
			hduout[0].header['COMMENT'] = '    {}*{} minus {}*{}'.format('%.2f'%a1, fn1.split('/')[-1], '%.2f'%a2, fn2.split('/')[-1])
			hduout.writeto(fnout, overwrite=overwrite)


	def _band_has_smaller_psf(self, band, bandto, fracdiff_threshold=0.00):
		"""
		whether band has smaller psf than bandto by a margin of "diffpsf_threshold" in arcsec. 
		get psf size from hsc_xid if the survey = 'hsc'. 

		Params
		------
		band (str)
		bandto (str)
		fracdiff_threshold=0.00:
			the threshold of psf difference in arcsec. If the diff is larger than return True. 

		Return
		------
		answer (bool)
		"""

		fp_stp = self.get_fp_stamp(band=band)
		fp_stpto = self.get_fp_stamp(band=bandto)

		return self._stamp_has_smaller_psf(fp_stp, fp_stpto, fracdiff_threshold=fracdiff_threshold)


	def _stamp_has_smaller_psf(self, fp_stp, fp_stpto, fracdiff_threshold=0.00, mode='moffat'):
		"""
		whether fp_stp has smaller psf than fp_stpto by a margin of "diffpsf_threshold" in arcsec. 
		get psf size from hsc_xid if the survey = 'hsc'. 

		Params
		------
		fp_stp (str):
			file path to stamp
		fp_stpto (str)
			file path to stamp to compare to 
		fracdiff_threshold=0.00:
			the threshold of psf difference in arcsec. If the diff is larger than return True. 

		Return
		------
		answer (bool)
		"""

		fp_psf = self.get_fp_psf(fp_stp)
		fp_psfto = self.get_fp_psf(fp_stpto)

		result = matchpsf.has_smaller_psf_fits(fp_psf, fp_psfto, mode=mode, fracdiff_threshold=fracdiff_threshold)

		return result


	def _get_psf_fwhm_pix(self, imgtag, mode='moffat'):
		""" return the moffat fwhm of the psf of the imgtag in pix"""
		fp_psf = self.get_fp_psf(self.get_fp_stamp_img(imgtag=imgtag))
		try:
			psf = fits.getdata(fp_psf)
			fwhm_pix = matchpsf.calc_psf_fwhm(psf, mode=mode)
		except:
			fwhm_pix = np.nan
		return fwhm_pix


	def _get_psf_fwhm_arcs(self, imgtag, mode='moffat'):
		""" return the moffat fwhm of the psf of the imgtag in arcsec"""
		fwhm_pix = self._get_psf_fwhm_pix(imgtag=imgtag, mode=mode)
		fwhm_arcs = fwhm_pix * np.array(self.pixsize.to(u.arcsec))
		return fwhm_arcs


	def _get_psf_moffat_params(self, imgtag):
		""" return gamma, alpha """
		fp_psf = self.get_fp_psf(self.get_fp_stamp_img(imgtag=imgtag))
		# try:
		psf = fits.getdata(fp_psf)
		model = matchpsf.fit_moffat(psf)
		gamma, alpha = model.gamma, model.alpha
		# except:
			# gamma, alpha = np.nan, np.nan
		return gamma, alpha


	def _copy_psf(self, fn_in, fn_out):
		""" given stamp paths copy the psf from in to out """
		fn_psf_in = self.get_fp_psf(fn_in)
		fn_psf_out = self.get_fp_psf(fn_out)
		copyfile(fn_psf_in, fn_psf_out)


