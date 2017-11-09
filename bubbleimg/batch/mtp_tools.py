"""
tools to make multiprocessing map work on class methods
https://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods
"""

import functools
import copyreg
import types

class partialmethod(functools.partial):
	def __get__(self, instance, owner):
		if instance is None:
			return self
		return functools.partial(self.func, instance, *(self.args or ()), **(self.keywords or {}))


def _pickle_method(method):
	func_name = method.__func__.__name__
	obj = method.__self__
	cls = method.__self__.__class__
	return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
	for cls in cls.mro():
		try:
			func = cls.__dict__[func_name]
		except KeyError:
			pass
		else:
			break
	return func.__get__(obj, cls)


copyreg.pickle(types.MethodType, _pickle_method, _unpickle_method)
