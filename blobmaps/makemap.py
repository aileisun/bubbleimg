 # subtractimg.py
# ALS 2015/08/10

"""
PURPOSE: To produce continuum subtracted broadband images and [OIII] intensity maps. It uses fromspec functions to determine the raito.  

"""
import os
import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as const


import fromspec
import imagedisp_util

from .. import filters

# import sys
# sys.path.append('../')
# import filters
# path_filters = os.path.dirname(filters.__file__)+'/'

nanomaggy = u.def_unit('nanomaggy', 3.631e-6*u.Jy)
u.add_enabled_units([nanomaggy])
u.nanomaggy=nanomaggy

def objw_makeall(obj, bandline='r',bandconti='z',update=False):
	"""
	PURPOSE: make all the image subtraction products. Including: 
		spec.pdf
		stamp-r-uconti.fits, stamp-r-gconti.fits, stamp-r-iconti.fits, stamp-r-zconti.fits 
		stamp-lOIII5008_F.fits
		stamp-lOIII5008_I.fits
	"""
	# check if end product exists
	files=[obj.dir_obj+'stamp-lOIII5008_F.fits',
			obj.dir_obj+'stamp-lOIII5008_I.fits',
			obj.dir_obj+'spec.fits',]
	
	isfiles=np.all([os.path.isfile(file) for file in files])

	if (not isfiles) or update:
		fromspec.plotSpecwFilters(obj)
		# # to compute all possible line minus conti
		# for band in ['u','g','r','i','z']:
		# 	if band != bandline:
		# 		objw_stamp_band_minusconti(obj, band1=bandline, band2=band, savefits=True)
		objw_stamp_band_minusconti(obj, band1=bandline, band2=bandconti, savefits=True)
		objw_stamp_lOIII5008_F(obj,bandline=bandline,bandconti=bandconti)
		objw_stamp_F_to_I(obj, mapname='lOIII5008', toplot=True)
		# objw_stamp_lOIII5008_I(obj)



def objw_makecontiscale(obj, bandline='r', bandconti='z', update=False):
	"""
	automatically make 'stamp-conti-onOIIIscale_F.fits' for obj
	"""
	files=[obj.dir_obj+'stamp-conti-onOIIIscale_F.fits',
			obj.dir_obj+'stamp-conti-onOIIIscale_I.fits']
	
	isfiles=np.all([os.path.isfile(file) for file in files])

	if (not isfiles) or update:
		objw_stamp_contionOIIIscale_F(obj, bandline=bandline,bandconti=bandconti)
		objw_stamp_F_to_I(obj, mapname='conti-onOIIIscale', toplot=False)


def objw_stamp_band_minusconti(obj, band1='r', band2='z', savefits=True):
	"""
	PURPOSE: subtract continuum in band1 using that scaled from band2

			 band1_contsub = band1 - contiratio12*band2

	PAREMETERS: 		 
			obj
			band1='r'
			band2='z'
			savefits=True

	RETURN: img (2d array)

	WRITE OUTPUT: obj.dir_obj+'stamp-'+band1+'-'+band2+'conti'+'.fits'
		          e.g., stamp-r-zconti.fits 
	"""
	# setup output filename
	filename=obj.dir_obj+'stamp-'+band1+'-'+band2+'conti'+'.fits'

	# get continuum ratio
	ratio12=fromspec.getObjBandContiRatio_dEdnu(obj, band1=band1, band2=band2)
	# subtract image
	img=subtractStampImg(obj, band1=band1, band2=band2, a1=1., a2=ratio12, savefits=True)

	# construct fits file
	header=fits.getheader(obj.dir_obj+'stamp-'+band1+'.fits')
	header['HISTORY']='Band continuum subtracted by ALS, see CONTIBAND, CONTIRTIO'
	header['CONTBAND']=band2
	header['CONTRTIO']=ratio12

	# write output
	prihdu = fits.PrimaryHDU(img, header=header)
	prihdu.writeto(filename, clobber=True)

	return img
	# d.set_np2arr(img)




def subtractStampImg(obj, band1='r', band2='z', a1=1., a2=1., savefits=True):
	"""
	PURPOSE: return subtracted stamp images (a1 * band1 - a2 * band2)

	INPUT: 
		obj (obj of class obsobj): an object containing information like RA, Dec, and SDSS PhotoObj. 
			Used attributes:  obj.dir_obj
		band1='r' (string)
		band2='z' (string)
		a1 (float) 			: normalization of image1
		a2 (float) 			: normalization of image2
		savefits (bool)

	READ INPUT: 
		obj.dir_obj+'stamp-'+band+'.fits'

	RETURN OUTPUT:
		image (2d array)

	WRITE OUTPUT: 
		obj.dir_obj+'stamp-'+band+'.fits'
	"""

	filename1=obj.dir_obj+'stamp-'+band1+'.fits'
	filename2=obj.dir_obj+'stamp-'+band2+'.fits'

	image1=fits.getdata(filename1)
	image2=fits.getdata(filename2)

	imagesub=a1*image1-a2*image2

	return imagesub

def objw_stamp_contionOIIIscale_F(obj, bandline='r',bandconti='z'):
	"""
	make a scaled continuum map from the bandconti, the scaling is such that 
	1 unit flux of continuum in this map has the same signal on the 
	'bandline' photometric image as the OIII line that has 1 unit flux. 

	maps are named 'stamp-conti-onOIIIscale_F.fits'
	"""

	# get z
	z=obj.sdss.z
	filein=obj.dir_obj+'stamp-'+bandconti+'.fits'
	fileout=obj.dir_obj+'stamp-conti-onOIIIscale_F.fits'

	# read in data
	hdu=fits.open(filein)
	img_bc=hdu[0].data*u.Unit(hdu[0].header['BUNIT'])

	# scale to the continuum level at the line band
	ratio12=fromspec.getObjBandContiRatio_dEdnu(obj, band1=bandline, band2=bandconti)
	img_bl=img_bc*ratio12

	# scale to the OIII5007 line intensity
	intR=filters.filtertools.intRoverldl(band=bandline)
	# get line ratios and wavelengths
	lOIII5008=5008.24*u.AA
	ROIII5007=fromspec.R_lambda(lOIII5008*(1.+z),band=bandline)

	img_scale=img_bl*intR*const.c/(ROIII5007*lOIII5008)
	contimap=img_scale.to(u.erg * u.s**-1 * u.cm**-2)

	# construct fits file
	hdu[0].data=contimap.value
	hdu[0].header['BUNIT']=(u.erg * u.s**-1 * u.cm**-2).to_string()
	hdu[0].header['SCALINE']='[OIII] 5007'
	hdu[0].header['SCABAND']=bandline
	hdu[0].header['FRMBAND']=bandconti
	hdu[0].header['HISTORY']='Converted from '+bandconti+' band image continuum level at the line intensity unit by ALS'

	# write fits file
	hdu.writeto(fileout, clobber=True)


def objw_stamp_lOIII5008_F(obj, bandline='r',bandconti='z'):
	"""
	PURPOSE: infer OIIIl5007 flux map [erg cm-2 s-1] from r-z subtracted images. 
	         and write file as obj.dir_obj+'stamp-lOIII5008_F.fits'

	DESCRIPTION: 

		The flux density measured in a band F_nu (like pixel values in images [nanomaggie ~ erg/cm2/s/Hz])
		are related to the spectrum f_nu(nu) and filter response funciton R(nu) by: 

			F_nu = integral{f_nu(nu)*R(nu)*dnu/nu}  / integral{R(nu)*dnu/nu}

			(Notice it's integrated over dnu/nu as CCD is counting photon, not energy. )

		If there are only two delta lines in the spectrum l1, l2, with fluxes f1, f2 we get:

			f2 = (F_nu * c * integral{R(lambda) * dlambda/lambda}) / (f1/f2 * R(l1)*l1 + R(l2)*l2)

		Now we have three lines Hb, [OIII]4959 [OIII]5007 and we want to get the flux of [OIII]5007

	WRITE OUTPUT: 
			obj.dir_obj+'stamp-lOIII5008_F.fits'
	"""
	from astropy.table import Table
	# setup
	filein=obj.dir_obj+'stamp-'+bandline+'-'+bandconti+'conti.fits'
	fileout=obj.dir_obj+'stamp-lOIII5008_F.fits'
	z=obj.sdss.z

	# sanity check - z range
	# the redshift must be within valid range where both [OIII] lines lie within r band (>60% of peak throughput)
	zranges=filters.filtertools.accessFile('OIIIredshiftrange0.6.txt')
	# zranges=Table.read('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/display/sdssdisp/test_filterrange/OIIIredshiftrange0.6.txt',format='ascii')
	z1,z2=zranges[zranges['band']==bandline]['z1','z2'][0]
	if z < z1 or z > z2: print ("WARNING: Invalid redshift. OIII falls outside filter 60 percent transmission range. ")

	else: 
		# read image of r-zconti with it's unit [nanomaggie]
		hdu=fits.open(filein)
		img=hdu[0].data*u.Unit(hdu[0].header['BUNIT'])

		# calculate integral{R(l) dl/l}
		intR=filters.filtertools.intRoverldl(band=bandline)
		# get line ratios and wavelengths
		lHb, lOIII4959, lOIII5008=4862.68*u.AA, 4959.295*u.AA, 5008.24*u.AA
		rHb_OIII5007= fromspec.getSDSSspeclineratio(obj,line1='H_beta',line2='[O_III] 5007')
		rOIII4959_OIII5007=1./2.98 # Theoretical calculation see Storey & Zeippen 00 2000MNRAS.312..813S
		# get filter response function at the wavelengths of the lines
		RHb,ROIII4959,ROIII5007=fromspec.R_lambda(lHb*(1.+z),band=bandline),fromspec.R_lambda(lOIII4959*(1.+z),band=bandline),fromspec.R_lambda(lOIII5008*(1.+z),band=bandline)

		# calculate line map
		numerator = img*intR*const.c	
		denominator = rHb_OIII5007*RHb*lHb*(1.+z)+rOIII4959_OIII5007*ROIII4959*lOIII4959*(1.+z)+ROIII5007*lOIII5008*(1.+z)
		linemap=(numerator/denominator).to(u.erg * u.s**-1 * u.cm**-2)

		# construct fits file
		hdu[0].data=linemap.value
		hdu[0].header['BUNIT']=(u.erg * u.s**-1 * u.cm**-2).to_string()
		hdu[0].header['LINE']='[OIII] 5007'
		hdu[0].header['HISTORY']='Converted from '+bandline+' band image to OIIIl5007 flux map by ALS'

		# write fits file
		hdu.writeto(fileout, clobber=True)
	


def objw_stamp_F_to_I(obj, mapname='lOIII5008', toplot=True):

	"""
	PURPOSE: infer OIIIl5007 intensity map [erg cm-2 s-1 arcsec-2] from flux maps stamp_F_lOIII5008 [erg cm-2 s-1]. 
			 see objw_stamp_F_lOIII5008 for more details.

	PAREMETERS: obj

	WRITE OUTPUT: 
			obj.dir_obj+'stamp-lOIII5008_I.fits'
	"""
	# setup
	filein=obj.dir_obj+'stamp-'+mapname+'_F.fits'
	fileout=obj.dir_obj+'stamp-'+mapname+'_I.fits'
	pixsize=0.396 *u.arcsec # SDSS pixsize [arcsec]

	# operation
	if os.path.isfile(filein):
		hdu=fits.open(filein)
		img_I=hdu[0].data*u.Unit(hdu[0].header['BUNIT'])/pixsize**2

		# modifying hdu
		hdu[0].data=img_I.value
		hdu[0].header['BUNIT']=img_I.unit.to_string()
		hdu.writeto(fileout, clobber=True)
		if toplot:
			imagedisp_util.fits_to_image(fileout,formats='pdf',toshow=False)
			imagedisp_util.fits_to_image(fileout,formats='png',toshow=False)
	else:
		print 'skipping making stamp_'+mapname+'_I.fits as input file stamp_'+mapname+'_F.fits does not exist'




# def objw_stamp_lOIII5008_I(obj):
# 	"""
# 	PURPOSE: infer OIIIl5007 intensity map [erg cm-2 s-1 arcsec-2] from flux maps stamp_F_lOIII5008 [erg cm-2 s-1]. 
# 			 see objw_stamp_F_lOIII5008 for more details.

# 	PAREMETERS: obj

# 	WRITE OUTPUT: 
# 			obj.dir_obj+'stamp-lOIII5008_I.fits'
# 	"""
# 	# setup
# 	filein=obj.dir_obj+'stamp-lOIII5008_F.fits'
# 	fileout=obj.dir_obj+'stamp-lOIII5008_I.fits'
# 	pixsize=0.396 *u.arcsec # SDSS pixsize [arcsec]

# 	# operation
# 	if os.path.isfile(filein):
# 		hdu=fits.open(filein)
# 		img_I=hdu[0].data*u.Unit(hdu[0].header['BUNIT'])/pixsize**2

# 		# modifying hdu
# 		hdu[0].data=img_I.value
# 		hdu[0].header['BUNIT']=img_I.unit.to_string()
# 		hdu.writeto(fileout, clobber=True)
# 		imagedisp_util.fits_to_image(fileout,formats='pdf',toshow=False)
# 		imagedisp_util.fits_to_image(fileout,formats='png',toshow=False)
# 	else:
# 		print 'skipping making stamp_lOIII5008_I.fits as input file stamp_lOIII5008_F.fits does not exist'


# def objw_lOIII5008_I_png(obj):
# 	"""
# 	to convert 'stamp-lOIII5008_I.fits' to 'stamp-lOIII5008_I.png' for obj
# 	"""
# 	fileout=obj.dir_obj+'stamp-lOIII5008_I.fits'
# 	imagedisp_util.fits_to_image(fileout,formats='png')

