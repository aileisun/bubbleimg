# fixture_built_batch.py
import os
import pytest
import shutil

from ..hscbatch import hscBatch

from setpaths import *


@pytest.fixture
def batch_good():
	if os.path.isdir(dir_batch):
		shutil.rmtree(dir_batch) 
	shutil.copytree(dir_verif, dir_batch)
	return hscBatch(dir_batch=dir_batch)


@pytest.fixture
def batch_wexcept():
	if os.path.isdir(dir_batch_wexcept):
		shutil.rmtree(dir_batch_wexcept) 
	shutil.copytree(dir_verif_wexcept, dir_batch_wexcept)
	return hscBatch(dir_batch=dir_batch_wexcept)


@pytest.fixture
def batch_onlyexcept():
	if os.path.isdir(dir_batch_onlyexcept):
		shutil.rmtree(dir_batch_onlyexcept) 
	shutil.copytree(dir_verif_onlyexcept, dir_batch_onlyexcept)
	return hscBatch(dir_batch=dir_batch_onlyexcept)


@pytest.fixture
def batch_confus():
	if os.path.isdir(dir_batch_confus):
		shutil.rmtree(dir_batch_confus) 
	shutil.copytree(dir_verif_confus, dir_batch_confus)
	return hscBatch(dir_batch=dir_batch_confus)

@pytest.fixture
def batch_hscphotoobj_incomplete():
	if os.path.isdir(dir_batch_hscphoto):
		shutil.rmtree(dir_batch_hscphoto) 
	shutil.copytree(dir_verif_hscphoto, dir_batch_hscphoto)
	return hscBatch(dir_batch=dir_batch_hscphoto, survey=survey)


@pytest.fixture
def batch_multiline():
	if os.path.isdir(dir_batch_multiline):
		shutil.rmtree(dir_batch_multiline) 
	shutil.copytree(dir_verif_multiline, dir_batch_multiline)
	return hscBatch(dir_batch=dir_batch_multiline, survey=survey, obj_naming_sys='sdss_precise')

