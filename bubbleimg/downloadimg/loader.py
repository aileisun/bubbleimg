# loader.py 
# ALS 2017/05/02

import numpy as np
import astropy.table as at
import astropy.units as u
import os 
import abc
from catalogue.catalogue_util import getSDSSName_fromRADEC

from ..filters import surveysetup
from ..class_obsobj import obsobj


class imgLoader(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, **kwargs):
		"""
		imgLoader

		Params
		----------
		/either
			obj (object of class obsobj): with attributes ra, dec, dir_obj
		/or  
			ra (float)
			dec (float)
			/either
				dir_obj (string)
			/or 
				dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'
				
		img_width (angle quantity or int) = 20 * u.arcsec:
			:if int then unit is assumed to be pix. 

		img_height (angle quantity or int) = 20 * u.arcsec:

		to_make_obj_sdss (bool) = False
			:if true, then make self.obj.sdss.xid etc and load xid.csv, photoboj.csv

		user = '': account user name for access
		password = '': account password for access

		Attributes
		----------
		ra (float)
		dec (float)
		dir_obj (string)
		img_width (angle quantity)
		img_height (angle quantity)
		img_width_pix (quantity of unit u.pix): floor integer of img_width in pixels
		img_height_pix (quantity of unit u.pix): floor integer of img_width in pixels

		obj (optional) may have attributes sdss.xid etc 

		"""
		#===== unparse input
		if 'obj' in kwargs: 
			self.obj = kwargs.pop('obj')
			self.ra = self.obj.ra
			self.dec = self.obj.dec
			self.dir_obj = self.obj.dir_obj
			
		else: 
			self.ra = kwargs.pop('ra', np.nan)
			self.dec = kwargs.pop('dec', np.nan)
			if 'dir_obj' in kwargs:
				self.dir_obj = kwargs.pop('dir_obj', '')
			elif 'dir_parent' in kwargs:
				sdssname = getSDSSName_fromRADEC(self.ra, self.dec)
				dir_parent = kwargs.pop('dir_parent', '')
				self.dir_obj = dir_parent+sdssname+'/'
			else:
				raise TypeError('dir_obj or dir_parent not specified')

		self.img_width = kwargs.pop('img_width', 20*u.arcsec)
		self.img_height = kwargs.pop('img_height', 20*u.arcsec)
		self.user = kwargs.pop('user', "")
		self.password = kwargs.pop('password', "")

		to_make_obj_sdss = kwargs.pop('to_make_obj_sdss', False)

		if kwargs:
			raise TypeError('Unexpected **kwargs: %r' % kwargs)

		if np.isnan(self.ra) or np.isnan(self.dec) or (self.dir_obj==''):
			raise TypeError('ra or dec or dir_obj not specified')
		#===== 

		# add on attributes
		if to_make_obj_sdss:
			self._make_obj_sdss_xid()

		self._attach_img_widthheight_unit()
		self.survey = 'to be overwritten'
		self.pixsize = -1


	@abc.abstractmethod
	def make_stamps(self, **kwargs):
		raise NotImplementedError("Subclass must implement abstract method")


	def get_stamp_filepath(self, band):
		return self.dir_obj+'stamp-{0}.fits'.format(band)


	def get_stamp_filename(self, band):
		return 'stamp-{0}.fits'.format(band)


	def _imgLoader__make_stamp_core(self, func_download_stamp, **kwargs):
		"""
		make a stamp of self and of band
		call _download_stamp and takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

		Params
		----------
		func_download_stamp (function)
		band (string) = 'r'
		overwrite (boolean) = False
		**kwargs to be entered into func_download_stamp()

		Return
		----------
		status: True if downloaded or skipped, False if download fails
		"""
		band = kwargs.pop('band', 'r')
		overwrite = kwargs.pop('overwrite', False)

		# setting
		filename = self.get_stamp_filename(band)
		file = self.get_stamp_filepath(band)

		if not os.path.isdir(self.dir_obj):
		    os.makedirs(self.dir_obj)

		if (not os.path.isfile(file)) or overwrite:
			print "download_stamp() ".format(filename)
			return func_download_stamp(band=band, **kwargs) # if failed then return False
		else:
			print "skip download_stamp() as file {0} exists".format(filename)
			return True


	def _imgLoader__make_stamps_core(self, func_download_stamp, **kwargs):
		"""
		make all stamp images of all survey bands for obj self. 

		Params
		----------
		func_download_stamp (function)
		overwrite (boolean) = False
		**kwargs to be entered into func_download_stamp()

		Return
		----------
		status: True if all downloaded or skipped, False if any of the downloads fails
		"""
		bands = surveysetup.surveybands[self.survey]

		statuss = np.ndarray(5, dtype=bool)
		for i, band in enumerate(bands): 
			statuss[i] = self._imgLoader__make_stamp_core(func_download_stamp=func_download_stamp, band=band, **kwargs)

		return all(statuss)


	def _attach_img_widthheight_unit(self):
		""" 
		make sure that self.img_width and self.height has angular unit 
		If they are a int/whole number then attach unit u.pix. otherwise raise error. 
		"""
		if not hasattr(self.img_width, 'unit'):
			if type(self.img_width) is int:
				self.img_width *= u.pix
				print "Assuming img_width unit is pix"
			elif type(self.img_width) is float and (self.img_width).is_integer():
				self.img_width *= u.pix
				print "Assuming img_width unit is pix"
			else: 
				raise ValueError("Param img_width has to be angular quantity or a whole number for pix")

		if not hasattr(self.img_height, 'unit'):
			if type(self.img_height) is int:
				self.img_height *= u.pix
				print "Assuming img_height unit is pix"
			elif type(self.img_height) is float and (self.img_height).is_integer():
				self.img_height *= u.pix
				print "Assuming img_height unit is pix"
			else: 
				raise ValueError("Param img_height has to be angular quantity or a whole number for pix")


	def _add_attr_img_width_pix_arcsec(self):
		""" creating attributes 
		self.img_width_pix
		self.img_height_pix
		self.img_height_arcsec
		self.img_width_arcsec
		which are self.img_width and img_height but changed to indicated units """
		survey_pixelscale = u.pixel_scale(self.pixsize/u.pixel)

		if hasattr(self.img_width, 'unit'): 
			self.img_width_pix = np.floor(self.img_width.to(u.pix, survey_pixelscale))
			self.img_width_arcsec = self.img_width.to(u.arcsec, survey_pixelscale)
		else: 
			raise ValueError('self.img_width has no angular units')

		if hasattr(self.img_height, 'unit'): 
			self.img_height_pix = np.floor(self.img_height.to(u.pix, survey_pixelscale))
			self.img_height_arcsec = self.img_height.to(u.arcsec, survey_pixelscale)
		else: 
			raise ValueError('self.img_height has no angular units')


	def _make_obj_sdss_xid(self):
		"""
		make sure that self have attributes self.obj.sdss.xid by calling _add_obj_sdss(update=True) when nessisary, and save files xid.csv, photoobj.csv. 
		"""
		try:
			self.obj.sdss.xid
		except(AttributeError):
			self._add_obj_sdss(update=True)			


	def _add_obj_sdss(self, update=True):
		""" 
		add self.obj with attributes ra, dec, sdss, sdss.xid, sdss.photoobj, etc., 
		download xid.csv, photoobj.csv into directory dir_obj 
		"""
		if (not hasattr(self, 'obj')) or update:
			dir_parent = '/'.join(self.dir_obj.split('/')[0:-2])+'/'
			tab = at.Table([[self.ra], [self.dec]], names=['ra', 'dec'])

			obj = obsobj(tab, catalog='sdss', dir_parent=dir_parent, towriteID=True)

			self.obj = obj

			if self.obj.dir_obj != self.dir_obj: 
				raise ValueError('self.dir_obj inconsistent with SDSS naming convension')

