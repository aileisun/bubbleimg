# linefrac.py

""" 
calculate the fractional contributions of lines in a given band, if all the flux in the band is contributed by lines 
"""
import numpy as np
import astropy.constants as const
import astropy.units as u

def fline_over_fnuband(fl, wl, Tl, fs, ws, Ts):
	"""
	calcualte the ratio flux_line / Fnu_band,

		{c / (wl * Tl)}  *  {( fl * wl * Tl) / sum[ fs * ws * Ts]  }

	where fl, wl, Tl are the flux, wavelength, and normalized filter transmission func of the target line. 
	where fs, ws, Ts are the flux, wavelength, and normalized filter transmission func of all the lines in band. 

	Params
	------
	fl (float): 
		flux of line
	wl (float)
		wavelength of line
	Tl (float):
		normalized filter transmission function at line

	\\ followings are the same but array of all the lines in band
	fs (array of float)
	ws (array of float)
	Ts (array of float)

	Return
	------
	ratio (quantity of unit Hz)
	"""

	# sanity check
	check_no_unit(fl)
	check_no_unit(fs)
	check_no_unit(Tl)
	check_no_unit(Ts)
	wl = check_either_no_unit_or_specificunit(wl, unit=u.AA)
	ws = check_either_no_unit_or_specificunit(ws, unit=u.AA)

	# calculation
	firstterm = const.c / (Tl * wl * u.AA)
	print("fs", fs)
	print("ws", ws)
	print("Ts", Ts)
	secondterm = FWT(fl, wl, Tl) / sumFWT(fs, ws, Ts)

	return firstterm * secondterm


def check_no_unit(quan):
	try:
		quan.unit
	except:
		return quan
	else: 
		raise Exception("Expecting no unit")


def check_either_no_unit_or_specificunit(quan, unit=u.AA):
	try:
		quan.unit
	except:
		return quan
	else: 
		if quan.unit != unit: 
			raise Exception("wavelength unit incorrect")
		else:
			return np.array(quan)


def FWT(f, w, T):
	return f*w*T


def sumFWT(fs, ws, Ts):
	return np.sum( fs * ws * Ts)