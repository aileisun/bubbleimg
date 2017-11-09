# plainobj.py
# ALS 2017/05/11

"""
define class plainObj, which has only ra, dec, and dir_obj as attributes.  
"""

import os

from .objnaming import get_obj_name


class plainObj(object):
	def __init__(self, ra=None, dec=None, obj_naming_sys='sdss', checkname=False, **kwargs):
		"""
		plainObj
		a object with only attributes ra, dec, name, and dir_obj. 
		If ra, dec is not provided then the name is taken from the directory of dir_obj. 
		it does not automatically create dir_obj when instantiated, one can call funciton make_dir_obj()


		Params
		----------
		ra = None (float)
		dec = None (float)
		obj_naming_sys = 'sdss' (string): 
			the naming system of object
		checkname = False (bool):
			whether to check if the directory name is consistent with ra, dec, given obj_naming_sys. 
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

		dir_parent (optional) (string)
		obj_naming_sys


		Note
		----
		It does not check if additional kwargs are inserted, as the following lines are not implemented. 
		if kwargs:
			raise TypeError('Unexpected **kwargs: %r' % kwargs)
		"""
		# if isinstance(ra, float):
		# 	self.ra = ra
		# else: 
		# 	raise TypeError('[plainobj] ra should be float')

		# if isinstance(dec, float):
		# 	self.dec = dec
		# else: 
		# 	raise TypeError('[plainobj] dec should be float')
		self.ra = ra
		self.dec = dec

		self.obj_naming_sys = obj_naming_sys

		obj_name = get_obj_name(self.ra, self.dec, obj_naming_sys=self.obj_naming_sys)

		if 'dir_obj' in kwargs:
			self.dir_obj = kwargs.pop('dir_obj', None)
			self.name = self.dir_obj.split('/')[-2]

			# sanity check: dir_obj naming consistent with ra, dec
			if checkname and (self.name[:4]=='SDSS' and self.name != obj_name):
				raise Exception('[plainobj] dir_obj SDSS name inconsistent with ra dec. {} != {}'.format
					(self.name, obj_name))
		elif 'dir_parent' in kwargs:
			dir_parent = kwargs.pop('dir_parent', None)
			self.name = obj_name
			self.dir_obj = dir_parent+self.name+'/'
			self.dir_parent = dir_parent
		else:
			raise Exception('[plainobj] dir_obj or dir_parent not specified')


		# sanity check: dir_obj is like a directory
		if self.dir_obj[-1] != '/':
			raise Exception('[plainobj] dir_obj not a directory path')



	def make_dir_obj(self):
		if not os.path.isdir(self.dir_obj):
			os.makedirs(self.dir_obj)