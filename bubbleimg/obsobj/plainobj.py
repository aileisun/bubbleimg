# plainobj.py
# ALS 2017/05/11

"""
define class plainObj, which has only ra, dec, and dir_obj as attributes.  
"""

import os

from objnaming import get_obj_name


class plainObj(object):
	def __init__(self, ra, dec, obj_naming_sys='sdss', **kwargs):
		"""
		plainObj
		a object with only attributes ra, dec, and dir_obj
		it does not automatically create dir_obj when instantiated, one can call funciton make_dir_obj()


		!!!! WARNING !!!!!
		no checking if unexpected kwargs is passed, in order to allow child class to have additional args
		if kwargs:
			raise TypeError('Unexpected **kwargs: %r' % kwargs)


		Params
		----------
		ra (float)
		dec (float)
		obj_naming_sys = 'sdss' (string): the naming system of object
		/either
			dir_obj (string)
		/or 
			dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'
				

		Attributes
		----------
		ra (float)
		dec (float)
		name (string)
		dir_obj (string)
		"""

		if isinstance(ra, float):
			self.ra = ra
		else: 
			raise TypeError('ra should be float')

		if isinstance(dec, float):
			self.dec = dec
		else: 
			raise TypeError('dec should be float')

		if 'dir_obj' in kwargs:
			self.dir_obj = kwargs.pop('dir_obj', None)
			self.name = self.dir_obj.split('/')[-2]

			# sanity check: dir_obj naming consistent with ra, dec
			if (self.name[:4]=='SDSS' and self.name != get_obj_name(self.ra, self.dec, obj_naming_sys='sdss')):
				raise Exception('dir_obj SDSS name inconsistent with ra dec')
		elif 'dir_parent' in kwargs:
			dir_parent = kwargs.pop('dir_parent', None)
			self.name = get_obj_name(self.ra, self.dec, obj_naming_sys=obj_naming_sys)
			self.dir_obj = dir_parent+self.name+'/'
			self.dir_parent = dir_parent
		else:
			raise Exception('dir_obj or dir_parent not specified')


		# sanity check: dir_obj is like a directory
		if self.dir_obj[-1] != '/':
			raise Exception('dir_obj not a directory path')



	def make_dir_obj(self):
		if not os.path.isdir(self.dir_obj):
			os.makedirs(self.dir_obj)