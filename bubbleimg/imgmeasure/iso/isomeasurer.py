# isomeasurer.py 
# ALS 2017/06/01
import os
import astropy.units as u
from astropy.io import fits
import numpy as np
# import matplotlib.pyplot as plt
import astropy.table as at
import pickle

from ..measurer import Measurer
import polytools
import plottools

class isoMeasurer(Measurer):

	def __init__(self, **kwargs):
		"""
		child of Measurer
		do isophotal measurements
		"""
		super(isoMeasurer, self).__init__(**kwargs)

		self.msrtype = 'iso'


	def get_fp_contours(self, imgtag='OIII5008_I', onlycenter=False):
		""" e.g., contours-{imgtag}.pkl 
			\or   contours_ctr-{imgtag}.pkl 
		"""
		if onlycenter:
			ctrtag = '_ctr'
		else:
			ctrtag = ''

		return self.dir_obj+'contours{ctrtag}-{imgtag}.pkl'.format(ctrtag=ctrtag, imgtag=imgtag)


	def make_measurements(self, imgtag='OIII5008_I', isocut=1.e-15*u.Unit('erg / (arcsec2 cm2 s)'), minarea=0, onlycenter=False, centerradius=2.*u.arcsec, overwrite=False, savecontours=False, plotmsr=False):
		"""
		make measurements on a map
			if imgtag='OIII5008_I' then measure 'stamp-OIII5008_I.fits'

		Params
		------
		self
		imgtag='OIII5008_I'
		overwrite = False (bool)

		isocut=1.e-15*u.Unit('erg / (arcsec2 cm2 s)'):
			isophote cut
		minarea=0:
			connected contour area (# pix) above the area is counted as part of the isophote measurement
		onlycenter=False:
			whether to consider only the center contours
		centerradius=2.*u.arcsec
		overwrite=False
		savecontours=False
		plotmsr=False

		Return
		------
		status (bool)

		Write Output 
		------------
		e.g., msr_iso-OIII5008.csv
		"""
		fn = self.get_fp_msr(imgtag=imgtag)

		if not os.path.isfile(fn) or overwrite:
			print("[isomeasurer] making measurement")

			fn_img = self.get_fp_stamp_img(imgtag=imgtag)
			# read in
			hdus = fits.open(fn_img)
			img = hdus[0].data * u.Unit(hdus[0].header['BUNIT'])
			xc, yc = self._get_xc_yc(img)

			# calc
			contours = self._get_contours_from_img(img=img, isocut=isocut, xc=xc, yc=yc, minarea=minarea, onlycenter=onlycenter, centerradius=centerradius)
			tab_msr = self._get_tab_measurements_from_contours(contours=contours, xc=xc, yc=yc)
			tab_params = self._get_tab_params(imgtag=imgtag, isocut=isocut, minarea=minarea, onlycenter=onlycenter, centerradius=centerradius)
			tabout = at.hstack([tab_params, tab_msr])

			# output
			tabout.write(fn, overwrite=overwrite)

			# optional output
			if savecontours:
				fn_contours = self.get_fp_contours(imgtag=imgtag, onlycenter=onlycenter)
				write_pickle(contours, fn_contours)

			if plotmsr:
				fn = self.get_fp_msrplot(imgtag=imgtag)
				fig, ax = plottools.plot_img(img, colorlabel=img.unit.to_string())
				plottools.overplot_contours(ax, contours)
				fig.savefig(fn)

		else:
			print("[isomeasurer] skip making measurement as files exist")

		return os.path.isfile(fn)


		
	def _get_tab_params(self, imgtag, isocut, minarea, onlycenter, centerradius):
		"""
		return a one row table of the measurement params 
		"""
		tab = at.Table([[imgtag], [str(isocut)], [minarea], [onlycenter], [str(centerradius)], ], names=['imgtag', 'isocut', 'minarea', 'onlycenter', 'centerradius', ])
		return tab


	def _get_tab_measurements_from_contours(self, contours, xc, yc):
		""" 
		calculate iso measurements from contours, return a table like: 
		"""

		tab = polytools.ShapeParamsTab_from_contours(contours, xc, yc)

		# unit conversion
		area_ars = tab['area_pix'][0]*(self.pixsize/u.arcsec)**2
		dmax_ars = self._pix_to_theta(tab['dmax_pix'][0], wunit=False)
		rmax_ars = self._pix_to_theta(tab['rmax_pix'][0], wunit=False)
		dper_ars = self._pix_to_theta(tab['dper_pix'][0], wunit=False)

		kpc_per_arcsec = np.array(self._get_kpc_proper_per_arcsec())

		area_kpc = area_ars * kpc_per_arcsec**2
		dmax_kpc = dmax_ars * kpc_per_arcsec
		rmax_kpc = rmax_ars * kpc_per_arcsec
		dper_kpc = dper_ars * kpc_per_arcsec

		tab_converted = at.Table(names=['area_kpc', 'dmax_kpc', 'rmax_kpc', 'dper_kpc', 'area_ars', 'dmax_ars', 'rmax_ars', 'dper_ars', ])
		tab_converted.add_row([area_kpc, dmax_kpc, rmax_kpc, dper_kpc, area_ars, dmax_ars, rmax_ars, dper_ars, ])

		tabout = at.hstack([tab_converted, tab])

		return tabout


	def _get_contours_from_img(self, img, isocut, xc, yc, minarea=0., onlycenter=False, centerradius=2.*u.arcsec):
		"""
		make contour at isocut of image as python pickle file (fn_contours)
		always overwrite

		Params
		------
		self
		img (array)
		isocut (float or quantity):
			has to be of the same type of unit as image
		minarea (float):
			minimum area (pix) to be considered as contour patch
		onlycenter (bool):
			whether to take only center patches as patches (they all have to pass minarea test as well)
		centerradius (angular quantity):
			if onlycenter = True, then what is the radius of the center area. only patches overlapping with that area will be considered. 

		"""

		# prep
		img_nparr = np.array((img/isocut).to(u.dimensionless_unscaled))

		# find contours -- satisfy minarea
		contours = polytools.find_largecontours(img=img_nparr, threshold=1., minarea=minarea)

		if onlycenter:  # select only those at the center
			centerradius_pix = self._theta_to_pix(centerradius)
			contours = polytools.select_center_contours(contours, xc, yc, radius=centerradius_pix)

		return contours



def read_pickle(fn):
    with open(fn, 'rb') as handle:
        result = pickle.load(handle)
    return result


def write_pickle(result, fn):
    with open(fn, 'wb') as handle:
        pickle.dump(result, handle)

