# plainobj.py
# ALS 2017/05/11

"""
define class plainObj, which has only ra, dec, and dir_obj as attributes.  
"""

import os

from catalogue.catalogue_util import getSDSSName_fromRADEC


class plainObj(object):
	def __init__(self, ra, dec, dirnamingsys='sdss', **kwargs):
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
		dirnamingsys = 'sdss' (string)
		/either
			dir_obj (string)
		/or 
			dir_parent (string): attr dir_obj is set to dir_parent+'SDSSJXXXX+XXXX/'
				

		Attributes
		----------
		ra (float)
		dec (float)
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
		elif 'dir_parent' in kwargs:
			dir_parent = kwargs.pop('dir_parent', None)
			sdssname = getSDSSName_fromRADEC(self.ra, self.dec)
			self.dir_obj = dir_parent+sdssname+'/'
			self.dir_parent = dir_parent
		else:
			raise Exception('dir_obj or dir_parent not specified')

		# sanity check: dir_obj is like a directory
		if self.dir_obj[-1] != '/':
			raise Exception('dir_obj not a directory path')

		# sanity check: dir_obj naming consistent with ra, dec
		dir_obj_name = self.dir_obj.split('/')[-2]
		if (dir_obj_name[:4]=='SDSS'):
			if (dir_obj_name!=getSDSSName_fromRADEC(self.ra, self.dec)): 
				raise Exception('dir_obj SDSS name inconsistent with ra dec')



	def make_dir_obj(self):
		if not os.path.isdir(self.dir_obj):
			os.makedirs(self.dir_obj)