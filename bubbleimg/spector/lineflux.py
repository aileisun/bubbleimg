# lineflux.py

"""
tool to estimate line flux and line flux error using trapz integral
"""
import numpy as np


def calc_line_flux(spec, ws, ivar, w0, w1, u_flux):
	""" calculate the flux and flux error of the line within the range w0 and w1 using trapz rule"""
	u_spec = spec.unit
	u_ws = ws.unit
	ivar = ivar.to(1./(u_spec**2))

	spec_uless = np.array(spec)
	ws_uless = np.array(ws)
	ivar_uless = np.array(ivar)

	if ivar.unit != (1./(spec.unit**2)):
		raise Exception("[spector] spec and ivar units inconsistent")

	# select region to integrate
	select_ws = (ws_uless > w0) & (ws_uless < w1)
	ws_sel = ws_uless[select_ws]
	spec_sel = spec_uless[select_ws]
	ivar_sel = ivar_uless[select_ws]
	var_sel = 1./ivar_sel

	# integrate
	f, fvar = trapz_var(x=ws_sel, y=spec_sel, yvar=var_sel)

	f = (f*u_spec*u_ws).to(u_flux)
	ferr = (np.sqrt(fvar)*u_spec*u_ws).to(u_flux)

	return f, ferr


def trapz(x, y):
	""" calculate the trapz integral of function y = f(x) """
	dx = np.diff(x)

	s = 0.

	for i in range(len(dx)):
		s += (y[i] + y[i+1])*0.5*dx[i]

	return s


def trapz_var(x, y, yvar):
	""" calculate the trapz integral of function y = f(x) and the variance of the integral"""
	dx = np.diff(x)

	s = 0.
	svar = 0.

	s = np.sum([(y[i] + y[i+1]) * 0.5*dx[i] for i in range(len(dx))])
	svar = np.sum([(yvar[i] + yvar[i+1]) * (0.5*dx[i])**2 for i in range(len(dx))])

	return s, svar
