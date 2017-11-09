
import pytest
from .. import getllambda
import imp
imp.reload(getllambda)

def test_getllambda():
	assert getllambda.getllambda('OIII5008') == 5008.240
	assert getllambda.getllambda('Ha') == 6564.61
	assert getllambda.getllambda('NII', lid=6585) == 6585.27
	assert getllambda.getllambda('OII', 3730) == 3729.88
	assert set(getllambda.getllambda('NII')) == set([5756.24, 6549.86, 6585.27])


def test_getllambda_err():
	with pytest.raises(Exception) as e:
		getllambda.getllambda('OIII5008', lid=4960)