# operator.py 
# ALS 2017/06/01

import abc

from .obsobj import obsObj


class Operator(object, metaclass=abc.ABCMeta):
	def __init__(self, **kwargs):
		"""
		objOperator

		Parent class for operators that operates on obsobj, e.g., imgLoader. 

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
				

		Attributes
		----------
		obj (instance of objObj)
		ra (float)
		dec (float)
		dir_obj (string)

		"""
		# unparse input		
		if 'obj' in kwargs: 
			self.obj = kwargs.pop('obj')
			# sanity check
			if 'dir_obj' in kwargs:
				if self.obj.dir_obj != kwargs.pop('dir_obj'):
					raise Exception("[operator] conflicting dir_obj entered")
		else: 
			self.obj = obsObj(**kwargs)

		self.ra = self.obj.ra
		self.dec = self.obj.dec
		self.dir_obj = self.obj.dir_obj

		# sanity check
		if self.dir_obj is None:
			raise TypeError('dir_obj not specified')
		