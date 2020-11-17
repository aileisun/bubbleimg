# isomeasurer.py 
# ALS 2017/06/01
import os
import astropy.units as u
from astropy.io import fits
import numpy as np
import astropy.table as at
import pickle
import scipy.ndimage as simg


from ..measurer import Measurer
from ... import tabtools
from . import polytools
from . import plottools

class isoMeasurer(Measurer):

	def __init__(self, **kwargs):
		"""
		child of Measurer
		do isophotal measurements
		"""
		super(isoMeasurer, self).__init__(**kwargs)

		self.msrtype = 'iso'


	def get_fp_contours(self, imgtag='OIII5008_I', onlycenter=False, suffix=''):
		""" e.g., msr_iso-OIII5008_I{suffix}_contours.pkl 
			\or   msr_iso-OIII5008_I{suffix}_contours-ctr.pkl 
		"""
		if onlycenter:
			ctrtag = '-ctr'
		else:
			ctrtag = ''

		fp_root = self.get_fp_msrtagroot(imgtag=imgtag, suffix=suffix)
		return fp_root+'_contours{ctrtag}.pkl'.format(ctrtag=ctrtag)


	def make_measurements(self, imgtag='OIII5008_I', isocut=3.e-15*u.Unit('erg / (arcsec2 cm2 s)'), minarea=5, onlycenter=True, centerradius=5.*u.arcsec, plotsuffix='', savecontours=False, plotmsr=False, msrsuffix='', overwrite=False, append=False):
		"""
		make measurements on a map and write to msr_iso.csv. 
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
		plotsuffix = '':
			plotsuffix label to be attach to the end of the plot or contour file names. 
		savecontours=False
		plotmsr=False
		msrsuffix=''
			plotsuffix label in the end of the measurement csv file: msr_iso_{msrsuffix}.csv.
		overwrite=False
		append=False

		Return
		------
		status (bool)

		Write Output 
		------------
		e.g., msr_iso.csv
		"""
		fn = self.get_fp_msr(msrsuffix=msrsuffix)

		condi = {'imgtag': imgtag, 'isocut': isocut, 'minarea': minarea, 'onlycenter': onlycenter, 'centerradius': centerradius}

		if append or overwrite or (not tabtools.fn_has_row(fn, condi)):
			print("[isomeasurer] making measurement")

			img = self.get_stamp_img(imgtag=imgtag, wunit=True)
			xc, yc = self._get_xc_yc(img)

			# calc
			if np.all(~np.isnan(img)): 
				contours = self._get_contours_from_img(img=img, isocut=isocut, xc=xc, yc=yc, minarea=minarea, onlycenter=onlycenter, centerradius=centerradius)
				tab_msr = self._get_tab_measurements_from_contours(contours=contours, xc=xc, yc=yc)
			else: 
				contours = []
				tab_msr = self._get_tab_measurements_nan()

			tab_params = self._get_tab_params(imgtag=imgtag, isocut=isocut, minarea=minarea, onlycenter=onlycenter, centerradius=centerradius)
			tabout = at.hstack([tab_params, tab_msr])

			# output
			tabtools.write_row(fn=fn, row=tabout, condi=condi, overwrite=overwrite, append=append)

			# optional output
			if savecontours:
				fn_contours = self.get_fp_contours(imgtag=imgtag, onlycenter=onlycenter, suffix=plotsuffix)
				write_pickle(contours, fn_contours, overwrite=overwrite)

			if plotmsr:
				fn_plot = self.get_fp_msrplot(imgtag=imgtag, suffix=plotsuffix)
				plottools.make_plot_img_w_contours(fn_plot=fn_plot, img=img, contours=contours)

		else:
			print("[isomeasurer] skip making measurement as files exist")

		return os.path.isfile(fn)


	def make_visualpanel(self, fn=None, compo_bands ='gri', imgtag='OIII5008_I', onlycenter=True, minarea=5, centerradius=5.*u.arcsec, tocolorbar=True, totitle=True, fontsize=12, overwrite=False):
		""" 
		make panel figure to visualize the composit and the iso measurements
		saved to e.g., 'msr_iso-OIII5008_I_panel.pdf'

		Params
		------
		fn = None: default: msr_iso_{imgtag}_panel.pdf
		compo_bands ='gri', imgtag='OIII5008_I', overwrite=False

		Return
		------
		status
		"""
		if fn is None:
			fn = self.get_fp_msrplot(imgtag=imgtag, suffix='_panel')
		else: 
			fn = self.dir_obj+fn

		if not os.path.isfile(fn) or overwrite:
			print("[isomeasurer] making visual panel")

			# get files ready 
			self.make_colorimg(bands=compo_bands, img_type='stamp', overwrite=False)

			# access data
			img_compo = simg.imread(self.dir_obj+'color_stamp-{}.png'.format(compo_bands))
			img_map = self.get_stamp_img(imgtag=imgtag, wunit=False)

			suffix = '_3e-15'
			isocut = 3.e-15*u.Unit('erg / (arcsec2 cm2 s)')
			fn_contours3 = self.get_fp_contours(imgtag=imgtag, onlycenter=onlycenter, suffix=suffix)
			if not os.path.isfile(fn_contours3):
				print("[isomeasurer] re-doing measurements to make contours required for visual panel plots")
				self.make_measurements(imgtag=imgtag, isocut=isocut, plotsuffix=suffix, minarea=minarea, onlycenter=onlycenter, centerradius=centerradius, overwrite=True, savecontours=True, plotmsr=False), 

			contours3 = read_pickle(fn_contours3)

			suffix = '_1e-15'
			isocut = 1.e-15*u.Unit('erg / (arcsec2 cm2 s)')
			fn_contours1 = self.get_fp_contours(imgtag=imgtag, onlycenter=onlycenter, suffix=suffix)
			if not os.path.isfile(fn_contours1):
				print("[isomeasurer] re-doing measurements to make contours required for visual panel plots")
				self.make_measurements(imgtag=imgtag, isocut=isocut, plotsuffix=suffix, minarea=minarea, onlycenter=onlycenter, centerradius=centerradius, overwrite=True, savecontours=True, plotmsr=False), 

			contours1 = read_pickle(fn_contours1)

			z = self.z
			pixsize = self.pixsize.to_value(u.arcsec)
			legend_suffix = ' at 3'
			name = self.obj.name[4:]

			title_compo = '${}~{}~{}~$'.format(compo_bands[0], compo_bands[1], compo_bands[2])+'$\mathrm{Composite}$'
			title_map = '$\mathrm{[OIII]\lambda 5007~Intensity}$'
			label_cbar = '$I~[10^{-15}~\mathrm{erg~s^{-1}~cm^{-2}~arcsec^{-2}}]$'

			plottools.make_iso_visual_panel(fn, img_compo, img_map, contours1, contours3, z, pixsize, legend_suffix, name, title_compo, title_map, label_cbar, tocolorbar=tocolorbar, totitle=totitle, fontsize=fontsize)

		else:
			print("[isomeasurer] skip making visual panel as files exist")

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


	def _get_tab_measurements_nan(self):
		""" 
		return a tab measurement just like _get_tab_measurements_from_contours() but with entries all nan. 
		"""
		names = ['area_kpc', 'dmax_kpc', 'rmax_kpc', 'dper_kpc', 'area_ars', 'dmax_ars', 'rmax_ars', 'dper_ars', 'area_pix', 'dmax_pix', 'rmax_pix', 'dper_pix', 'theta_dmax', 'theta_rmax', 'theta_dper', 'aspectr']

		tabout = at.Table(names=names)
		tabout.add_row([np.nan for i in range(len(names))])

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
			if onlycenter = True, then it sets the radius of the center area. only patches overlapping with that area will be considered. 

		"""

		# prep
		try: 
			img.unit
		except:
			img_nparr = img/isocut
		else:
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


def write_pickle(result, fn, overwrite=False):

	if not os.path.isfile(fn) or overwrite:
		with open(fn, 'wb') as handle:
			pickle.dump(result, handle)

