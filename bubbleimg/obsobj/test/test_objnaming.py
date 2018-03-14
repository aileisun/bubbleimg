# test_objnaming.py
# ALS 2017/07/20

import pytest

from .. import objnaming


ra = 150.0547735
dec = 12.7073027

def test_objnaming():

	name = objnaming.get_obj_name(ra, dec, obj_naming_sys='sdss')
	assert name == 'SDSSJ1000+1242'

	name = objnaming.get_obj_name(ra, dec, obj_naming_sys='sdss_precise')
	assert name == 'SDSSJ100013+124226'

	name = objnaming.get_obj_name(ra, dec, obj_naming_sys='j')
	assert name == 'J1000+1242'

	name = objnaming.get_obj_name(ra, dec, obj_naming_sys='j_precise')
	assert name == 'J100013+124226'

	name = objnaming.get_obj_name(ra=29.158592, dec=-4.0001336, obj_naming_sys='sdss_precise')
	assert name == 'SDSSJ015638-040000'

	name = objnaming.get_obj_name(ra=140.513341745004, dec=-0.745408930047981, obj_naming_sys='sdss_precise')
	assert name == 'SDSSJ092203-004443'

