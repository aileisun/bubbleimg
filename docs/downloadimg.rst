***************
Download Images
***************

The module downloadimg is for downloading and constructing stamp (cut-out) images. Currently only SDSS and HSC are supported. 


Download SDSS stamp images
==========================

- Specify Source
One ways is to download image given a ra, dec (decimal). 

	>>> from bubbleimg.downloadimg import SDSSimgLoader
	>>>
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_obj = './here/'
	>>> L = SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj)
	>>> status = L.make_stamps()

The fits file of the aligned five band SDSS images will be stored as './here/stamp-r.fits', etc. If status == True then the download is successful. 

Or one can also download image given an obsobj object that has the required ra, dec, dir_obj info. 

	>>> import astropy.table as at
	>>> from bubbleimg.downloadimg import SDSSimgLoader
	>>> from bubbleimg import obsobj

	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_parent = './'

	>>> tab = at.Table([[ra], [dec]], names=['ra', 'dec'])
	>>> obj = obsobj(tab, catalog='SDSS', dir_parent=dir_parent, towriteID=False)
	>>> L = SDSSimgLoader(obj=obj)
	>>> status = L.make_stamps()


- Image Size

The image size is specified with the img_width and img_height parameters. Numeric input will be taken as number of pixels. 

	>>> L = SDSSimgLoader(obj=obj, img_width=64, img_height=64)
	>>> status = L.make_stamps()

One can also specify the angular scale. 

	>>> import astropy.units as u
	>>> L = SDSSimgLoader(obj=obj, img_width=20*u.arcsec, img_height=20*u.arcsec)
	>>> status = L.make_stamps()


- Other Behavior

The default behavior of SDSSimgLoader is to not overwrite existing files, use 'r' band as reference band for alignment, and not to keep the frame images after stamps are done. These options can be changed with parameters ``band_rf`` and ``tokeepframe``. 

	>>> status = L.make_stamps(overwrite=False, band_rf='r', tokeepframe=False)


Download HSC stamp images
=========================

Downloading HSC images is largely the same as SDSS, except that you need to provide your STARs username and password. You can choose to download one band only. 

	>>> from bubbleimg.downloadimg import HSCimgLoader
	>>>
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_obj = './here/'
	>>> L = HSCimgLoader(ra=ra , dec=dec, dir_obj=dir_obj, user='XXXX', password='XXXXXXXX')
	>>> status = L.make_stamp(band='r')
	>>> status = L.make_stamps()
