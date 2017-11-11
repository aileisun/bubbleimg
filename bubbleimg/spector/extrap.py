# extrap.py
# 07/06/2017 ALS 

"""
tools to extrapolate spectrum. Currently works on continuum. 
"""

import numpy as np
from . import getconti

def extrapolate(ys1, xs1, xs2, polydeg=1, extbase_length=2000., epsilon=1.e-5):
	""" extrapolate ys1 to cover the entire range of xs2 """
	if not (min(xs1) <= min(xs2)):# ys1 is short at the low x end
		x_end = min(xs2)+epsilon

		ys1, xs1 = extrap_to_end(ys1, xs1, x_end = x_end, polydeg=polydeg, extbase_length=extbase_length)
	
	if not (max(xs1) >= max(xs2)):# ys1 is short at the high x end
		x_end = max(xs2)-epsilon

		ys1, xs1 = extrap_to_end(ys1, xs1, x_end = x_end, polydeg=polydeg, extbase_length=extbase_length)

	return ys1, xs1


def extrap_to_ends(ys, xs, x_end0, x_end1, polydeg=1, extbase_length=2000.):
	# extrapolte blue end
	ys_ext, xs_ext = extrap_to_end(ys=ys, xs=xs, x_end0=l0, polydeg=polydeg, extbase_length=extbase_length)
	# extrapolte red end
	ys_ext, xs_ext = extrap_to_end(ys=ys_ext, xs=xs_ext, x_end=l1, polydeg=polydeg, extbase_length=extbase_length)

	return ys_ext, xs_ext


def extrap_to_end(ys, xs, x_end, polydeg=1, extbase_length=2000.):
	""" 
	get ys_ext, xs_ext that is extrapolated ys, xs out to x_end.  

	Params
	------
	ys: ys to extrapolate
	xs: ys's coord system
	x_end: to extrapolate to

	Return
	------
	ys_ext, xs_ext	
	"""
	ys_uless = np.array(ys)
	xs_uless = np.array(xs)
	x_end = np.array(x_end)
	
	if not arr_enclose_x(xs_uless, x_end):
		print("[extrap] extrapolating")

		xs_ext, xs_add = extend_xs_to_end(xs_uless, x_end)

		ys_ext = _extrap_polyfit(ys_uless, xs_uless, x_end, xs_add, polydeg=polydeg, extbase_length=extbase_length)

		ys_ext = getconti.inherit_unit(ys_ext, ys)
		xs_ext = getconti.inherit_unit(xs_ext, xs)

		return ys_ext, xs_ext
	else:
		print("[extrap] skip extrapolating as xs contains x_end")
		return ys, xs


def _extrap_polyfit(ys_uless, xs_uless, x_end, xs_add, polydeg=1, extbase_length=2000.):
	""" 
	use polynomial fit to extrap the red end of the spectrum 

	Params
	------
	polydeg: deg of poly
	extbase_length: 
		the length on the side of array in x space to do poly fitting on in units of xs's unit (usually Anstrom). 


	"""
	ys_base, xs_base = get_extbase(ys_uless, xs_uless, x_end, extbase_length=extbase_length)
	param = np.polyfit(x=xs_base, y=ys_base, deg=polydeg)
	p = np.poly1d(param)
	ys_add = p(xs_add)
	x_start, trailing = get_trailing(xs_uless, x_end)
	if trailing:
		ys_ext = np.append(ys_uless, ys_add)
	else: 
		ys_ext = np.append(ys_add, ys_uless)

	return ys_ext


def arr_enclose_x(arr, x):

	return ((np.max(arr) >= x) and (np.min(arr) <= x))



def get_extbase(ys, xs, x_end, extbase_length):	
	""" get the red end of the spectrum to do estimation for extrapolation """
	x_start, trailing = get_trailing(xs, x_end)

	sel = [np.absolute(x-x_start) < extbase_length for x in xs]

	ys_base = ys[sel]
	xs_base = xs[sel]

	return ys_base, xs_base


def get_trailing(xs, x_end):
	""" 
	if x_end is larger than xs
	return x_start (the end of xs facing x_end) and trailing (bool) indicating whether x_end > xs
	complains if x_end is within xs
	"""
	if x_end > max(xs):
		x_start = max(xs)
		trailing = True

	elif x_end < min(xs):
		x_start = min(xs)
		trailing = False
	else:
		raise Exception("[extrap] x_end within xs")

	return x_start, trailing


def extend_xs_to_end(xs, x_end):
	""" 
	extend the array of xs to cover x_end 
	if xs is equal spacing in log/linear space, conserve the equal spacing. Otherwise, the additional array is equally spaced in linear space using averaged step size. 
	"""

	x_start, trailing = get_trailing(xs, x_end)

	if is_equalspace_log(xs):
		dlnx = get_dlnx(xs)
		xs_add = geomspace_stepsize(x_start, x_end, dlnx)[1:]

	else:
		dx = get_dx(xs)
		xs_add = linspace_stepsize(x_start, x_end, dx)[1:]

	xs_add = np.sort(xs_add)

	if trailing:
		xs_new = np.append(xs, xs_add)
	else: 
		xs_new = np.append(xs_add, xs)

	return xs_new, xs_add


def is_equalspace_log(xs, sigdigit=5):
	""" if array xs is equally spaced in log space"""
	return len(set(np.around(np.diff(np.log(xs)), sigdigit))) == 1


def get_dlnx(xs):
	return np.mean(np.diff(np.log(xs)))


def geomspace_stepsize(start, stop, dlnx):
	""" return equally log spaced array from start to stop (sometimes over) with step size dlnx """
	num = np.ceil(np.absolute(np.log(stop) - np.log(start))/dlnx) + 1

	if stop > start:
		stop_real = np.exp(np.log(start)+(num-1)*dlnx)
	else: 
		stop_real = np.exp(np.log(start)-(num-1)*dlnx)

	return np.geomspace(start, stop_real, num=num, endpoint=True)


def get_dx(xs):
	return np.mean(np.diff(xs))


def linspace_stepsize(start, stop, dx):
	""" return equally lin spaced array from start to stop (sometimes over) with step size dx """
	num = np.ceil(np.absolute(stop - start)/dx) + 1

	if stop > start:
		stop_real = start+(num-1)*dx
	else: 
		stop_real = start-(num-1)*dx

	
	return np.linspace(start, stop_real, num=num, endpoint=True)
