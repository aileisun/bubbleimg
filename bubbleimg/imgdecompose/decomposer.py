# decomposer.py 
# ALS 2017/06/01


from ..obsobj import Operator
from ..spectrum import spec

class Decomposer(Operator):

	def __init__(self, **kwargs):
		"""
		Decomposer, an image Operator

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

		survey = 'hsc'
				

		Attributes
		----------
		Operator Attributes:	
			obj (instance of objObj)
			ra (float)
			dec (float)
			dir_obj (string)

		survey = 'hsc' (string)
		"""
		#===== unparse input
		super(Decomposer, self).__init__(**kwargs)

		self.survey = kwargs.pop('survey', 'hsc')



	def make_spec_SEDs(self):
		"""
		make the continuum and line SEDs from spectrum 
		"""
		pass