
import numpy as np
import astropy.table as at
from astropy.io import fits
import os
import requests

from ... import external_links

def psField_to_psfs(dir_obj='./', photoobj=None, bands=['u', 'g', 'r', 'i', 'z']):
	""" 
	For a given dir_obj translate dir_obj/psField.fits to dir_obj/psf-(band).fits file using readatlas given the band, and take out the bias of 1000 and normalize it to one. 
	If photoboj not supplied, read dir_obj/PhotoObj.csv as photoobj table. 

	Params
	------
	dir_obj='./' (string)
	photoobj=None (astropy table): table with columns rowc_(band) calc_(band), e.g. obj.sdss.photobj
	bands=['u', 'g', 'r', 'i', 'z']
	"""
	infile = dir_obj+'psField.fits'

	if os.path.isfile(infile):
		for ib, band in enumerate(bands, 1):
			fn_raw = dir_obj+'psf-'+band+'_raw.fits'
			fn_psf = dir_obj+'psf-'+band+'.fits'
			psField_to_raw_psf(fn_psField=infile, fn_out=fn_raw, photoobj=photoobj, band=band, ib=ib)
			raw_psf_to_psf(fn_raw, fn_psf)

			os.remove(fn_raw)

	else: 
		raise ValueError('psField.fits does not exist')


def raw_psf_to_psf(fn_raw, fn_psf):
	"""
	translate raw psf file, which is biased by 1000 and not normalized to psf file which is not biased and is normalized. 
	"""
	hdus = fits.open(fn_raw)
	data = hdus[0].data
	data = data - 1000.
	data = data/np.sum(data)

	hdus[0].data = data

	if round(hdus[0].data[0,0], 5) != 0.:
		raise Exception("The edge of PSF is not close to zero")

	hdus.writeto(fn_psf, overwrite=True)


def psField_to_raw_psf(fn_psField, fn_out, photoobj, band, ib):
	"""
	translate psField to psf-(band), which has a bias of 1000 and is not normalized to 1. 
	"""
	if photoobj is None:
		photoobj = at.Table.read(dir_obj+'PhotoObj.csv',format='ascii.csv')	

	row = photoobj['rowc_'+band]
	col = photoobj['colc_'+band]

	# running command
	command = get_readpsf_command(fn_psField, ib, row, col, fn_out)
	os.system(command)


def get_readpsf_command(infile, ib, row, col, outfile):
	""" get the command for reading psf """
	codepath = external_links.dir_readatlas
	return codepath+'read_PSF '+infile+' '+'%d'%ib+' '+'%.4f'%row+' '+'%.4f'%col+' '+outfile


def download_psField(xid, dir_out='./', filename='psField.fits'):
	"""
	download the psField to directory dir_out as 'psField.fits'

	Params
	------
	xid: astropy table with columns ['run', 'rerun', 'camcol', 'field']]
	dir_out='./' (string)
	filename='psField.fits' (string)
	"""
	print("downloading psField file for psf")

	# sanity check
	if not isinstance(xid, at.Table) or len(xid)!=1:
		raise ValueError("xid is not a one row table as needed")

	filepath = dir_out+filename

	# constructing URL
	run, rerun, camcol, field = [xid[key][0] for key in ['run', 'rerun', 'camcol', 'field']]
	url = get_url_psField(run, rerun, camcol, field)

	# download
	rqst = requests.get(url)
	if rqst.status_code == 200:
		with open(filepath, 'wb') as out:
			for bits in rqst.iter_content():
				out.write(bits)
		return True
	else:  
		print("psField cannot be retrieved")
		return False


def get_url_psField(run, rerun, camcol, field):
	urlsdssbase='http://data.sdss3.org/sas/dr12/'
	urlpsField='boss/photo/redux/'+'%d'%rerun+'/'+'%d'%run+'/objcs/'+'%d'%camcol+'/psField-'+'%06d'%run+'-'+'%d'%camcol+'-'+'%04d'%field+'.fit'
	return urlsdssbase+urlpsField

