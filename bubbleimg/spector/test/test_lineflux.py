import numpy as np

from .. import lineflux


def test_trazperr_trapz():
	x = np.arange(10)
	y = x**2

	trapz = np.trapz(y=y, x=x)

	my_trapz = lineflux.trapz(x=x, y=y)

	assert np.absolute(my_trapz-trapz)/trapz < 1.e-5


def test_trazperr_trapz_var():
	x = np.array([0., 3.])
	y = np.array([0., 1.])
	yvar = np.array([1., 1.])

	integral, variance = lineflux.trapz_var(x=x, y=y, yvar=yvar)

	assert integral == 1.5
	assert variance == np.sum(yvar**2)*(1.5**2)
