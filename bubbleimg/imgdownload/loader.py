# loader.py 
# ALS 2017/05/02

import numpy as np
import astropy.table as at
import astropy.units as u
import os 
# import abc

from ..filters import surveysetup
from .. import obsobj
from ..external_links import file_humvi_compose


class imgLoader(obsobj.Operator):
	# __metaclass__ = abc.ABCMeta

	def __init__(self, **kwargs):
		"""
		imgLoader, child of Operator

		Params
		----------
		Operator params:
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

		user = '': account user name for access
		password = '': account password for access

		Attributes
		----------
		Operator Attributes:	
			obj (instance of objObj)
			ra (float)
			dec (float)
			dir_obj (string)
		status:
			whether loader is properly initiated, for example, successfully queried xid. If False, than make file functions will not be executed and returns False. 
		img_width (angle quantity)
		img_height (angle quantity)
		img_width_pix (quantity of unit u.pix): floor integer of img_width in pixels
		img_height_pix (quantity of unit u.pix): floor integer of img_width in pixels

		Attributes to be overwritten by subclass:
			survey (str)
			pixsize (astropy quantity in arcsec)
		"""
		#===== unparse input
		super(imgLoader, self).__init__(**kwargs)

		self.status = False # to be overwritten by subclass

		self.img_width = kwargs.pop('img_width', 20*u.arcsec)
		self.img_height = kwargs.pop('img_height', 20*u.arcsec)
		self._user = kwargs.pop('user', "")
		self._password = kwargs.pop('password', "")
		
		self._attach_img_widthheight_unit()
		self.survey = 'to be overwritten'
		self.pixsize = -1


	# @abc.abstractmethod
	def make_stamps(self, **kwargs):
		raise NotImplementedError("Subclass must implement abstract method")


	# @abc.abstractmethod
	def make_psfs(self, **kwargs):
		raise NotImplementedError("Subclass must implement abstract method")


	def get_fp_stamp(self, band):
		return self.dir_obj + self.get_fn_stamp(band)


	def get_fn_stamp(self, band):
		return 'stamp-{0}.fits'.format(band)


	def get_fp_psf(self, band):
		return self.dir_obj + self.get_fn_psf(band)


	def get_fn_psf(self, band):
		return 'psf-{0}.fits'.format(band)


	def _imgLoader__make_file_core(self, func_download_file, func_naming_file, band='r', overwrite=False, **kwargs):
		"""
		make a file of self and of band
		call _download_file and takes care of overwrite with argument 'overwrite'. Default: do not overwrite. 

		Params
		----------
		func_download_file (function), which takes (self, band) as argument
		func_naming_file (function), which takes (self, band) as argument
		band (string) = 'r'
		overwrite (boolean) = False
		**kwargs to be entered into func_download_file()

		Return
		----------
		status: True if downloaded or skipped, False if download fails
		"""

		# setting
		if self.status:
			filename = func_naming_file(band)
			filepath = self.dir_obj+filename

			if not os.path.isdir(self.dir_obj):
			    os.makedirs(self.dir_obj)

			if (not os.path.isfile(filepath)) or overwrite:
				print("[loader] making file {}".format(filename))
				# try: 
				return func_download_file(band=band, **kwargs) # if failed then return False
				# except:
					# return False
			else:
				print("[loader] skip make file as file {0} exists".format(filename))
				return True
		else: 
			print("[loader] make file failed as loader not properly initiated")
			return False


	def _imgLoader__make_files_core(self, func_download_file, func_naming_file, overwrite=False, **kwargs):
		"""
		make all file images of all survey bands for obj self. 

		Params
		----------
		func_download_file (function), which takes (self, band) as argument
		func_naming_file (function), which takes (self, band) as argument
		overwrite (boolean) = False
		**kwargs to be entered into func_download_file()

		Return
		----------
		status: True if all downloaded or skipped, False if any of the downloads fails
		"""
		bands = surveysetup.surveybands[self.survey]

		statuss = np.ndarray(len(bands), dtype=bool)
		for i, band in enumerate(bands): 
			statuss[i] = self._imgLoader__make_file_core(func_download_file=func_download_file, func_naming_file=func_naming_file, band=band, overwrite=overwrite, **kwargs)

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
		self.img_width_pix (int)
		self.img_height_pix (int)
		self.img_height_arcsec (angular quantity)
		self.img_width_arcsec (angular quantity)
		which are self.img_width and img_height but changed to indicated units """
		survey_pixelscale = u.pixel_scale(self.pixsize/u.pixel)

		if hasattr(self.img_width, 'unit'): 
			nwpix = (self.img_width.to(u.pix, survey_pixelscale)/u.pix).to(u.dimensionless_unscaled)
			self.img_width_pix = int(np.floor(nwpix))
			self.img_width_arcsec = self.img_width.to(u.arcsec, survey_pixelscale)
		else: 
			raise ValueError('self.img_width has no angular units')

		if hasattr(self.img_height, 'unit'): 
			nhpix = (self.img_height.to(u.pix, survey_pixelscale)/u.pix).to(u.dimensionless_unscaled)
			self.img_height_pix = int(np.floor(nhpix))
			self.img_height_arcsec = self.img_height.to(u.arcsec, survey_pixelscale)
		else: 
			raise ValueError('self.img_height has no angular units')


	def add_obj_sdss(self, update=False):
		""" 
		if self.obj does not exist create self.obj as an instance of obsObj
		if self.obj does not have sdss then call add_sdss(), which makes self.obj.sdss.xid, etc. 
		download xid.csv, photoobj.csv into directory dir_obj 

		Return
		------
		status: True if success, False if not
		"""
		if (not hasattr(self, 'obj')) or update:
			self.obj = obsobj.obsObj(ra=self.ra, dec=self.dec, dir_obj=self.dir_obj)

		if (not hasattr(self.obj, 'sdss')) or update:
			status = self.obj.add_sdss()
		else: 
			status = True

		# sanity check
		if self.obj.dir_obj != self.dir_obj: 
			raise ValueError('self.dir_obj inconsistent with SDSS naming convension')

		return status


	def add_obj_hsc(self, update=False, **kwargs):
		""" 
		if self.obj does not exist create self.obj as an instance of obsObj
		if self.obj does not have hsc then call add_hsc(), which makes self.obj.hsc.xid, etc. 
		download xid.csv, photoobj.csv into directory dir_obj 

		Return
		------
		status: True if success, False if not
		"""
		if (not hasattr(self, 'obj')) or update:
			self.obj = obsobj.obsObj(ra=self.ra, dec=self.dec, dir_obj=self.dir_obj)

		if (not hasattr(self.obj, 'hsc')) or update:
			status = self.obj.add_hsc(**kwargs)
		else:
			status = True

		# sanity check
		if self.obj.dir_obj != self.dir_obj: 
			raise ValueError('self.dir_obj inconsistent with SDSS naming convension')

		return status


	def plot_colorimg(self, bands ='riz', img_type='stamp', overwrite=False):
		"""
		make color composit image using external package HumVI. Example file name: 'color_stamp-riz.png'.

		Params
		------
		bands ='riz'
		img_type='stamp'
		overwrite=False

		Return
		------
		status (bool)
		"""
		fn = self.dir_obj+'color_{img_type}-{bands}.png'.format(bands=bands, img_type=img_type)

		fns_in = [self.dir_obj+img_type+'-'+band+'.fits' for band in bands[::-1]]

		if (not os.path.isfile(fn)) or overwrite:
			commandfiles = '{0} {1} {2} {3}'.format(fn, fns_in[0], fns_in[1], fns_in[2])
			commandHumVI = file_humvi_compose+' -s 1.0,1.1,1.0  -p 1.6,1.6  -o '+commandfiles

			os.system(commandHumVI)

		status = os.path.isfile(fn)

		return status

