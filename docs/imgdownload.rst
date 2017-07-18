***************
Download Images
***************

The module imgdownload is for downloading and constructing stamp (cut-out) images. Currently only SDSS and HSC are supported. 


Download SDSS stamp images
==========================

Specify Source
--------------

One ways is to provide a ``ra``, ``dec`` pair (decimal in degree) and a directory ``dir_obj``. 

	>>> from bubbleimg.imgdownload import SDSSimgLoader
	>>>
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_obj = './here/'
	>>> L = SDSSimgLoader(ra=ra , dec=dec, dir_obj=dir_obj)


Or one can also provide an obsobj object that has the required info. 

	>>> import astropy.table as at
	>>> from bubbleimg.imgdownload import SDSSimgLoader
	>>> from bubbleimg import obsobj
	>>>
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_obj = './here/'
	>>>
	>>> obj = obsobj(ra=ra, dec=dec, dir_obj=dir_obj)
	>>> L = SDSSimgLoader(obj=obj)


``SDSSimgLoader`` will automatically try to find a matching SDSS object close to the coordinate. One can check whether it succeed by:

	>>> L.status


Download stamps and psf
-----------------------
Once you have an ``SDSSimgLoader`` defined, you can download the files by.

	>>> status = L.make_stamps()

The fits file of the aligned five band SDSS images will be stored as './here/stamp-r.fits', etc. If status == True then the download is successful. 

Similarly, one can download psfs by:

	>>> status = L.make_psfs()


Data Release
------------

The default data release is set to DR13, one can specify the data release by, for example: 

	>>> L = SDSSimgLoader(obj=obj, data_release=7)


Image Size
----------

The image size is specified with the ``img_width`` and ``img_height`` parameters. Numeric input will be taken as number of pixels. 

	>>> L = SDSSimgLoader(obj=obj, img_width=64, img_height=64)
	>>> status = L.make_stamps()

One can also specify the angular scale. 

	>>> import astropy.units as u
	>>> L = SDSSimgLoader(obj=obj, img_width=20*u.arcsec, img_height=20*u.arcsec)
	>>> status = L.make_stamps()


Other Behaviors
---------------

The default behavior of ``SDSSimgLoader`` is to not overwrite existing files, use 'r' band as reference band for alignment, and not to keep the frame images after stamps are done. These options can be changed with parameters ``band_rf`` and ``tokeepframe``. 

	>>> status = L.make_stamps(overwrite=False, band_rf='r', tokeepframe=False)



Download HSC stamp images
=========================

STARs account
-------------

Downloading HSC images is largely the same as SDSS, except that you need to provide your STARs username and password. The easiest way is to set your STARS username and passwords as environmental variable. In the command line:

	>>> export HSC_SSP_CAS_USERNAME
	>>> read -s HSC_SSP_CAS_USERNAME
	>>> export HSC_SSP_CAS_PASSWORD
	>>> read -s HSC_SSP_CAS_PASSWORD

Otherwise you will be asked to enter it manually each time bubbleimg access HSC data. 


Specify Source
--------------

Just like SDSS, to initiate and HSC, you can use parameters (``ra``, ``dec``, ``dir_obj``) or an ``obsobj`` instance as parameter. For example if you have an obsobj

	>>> from bubbleimg import obsobj
	>>> 
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_obj = './here/'
	>>> obj = obsobj.obsObj(ra=ra, dec=dec, dir_obj=dir_obj)

then

	>>> from bubbleimg.imgdownload import HSCimgLoader
	>>> L = HSCimgLoader(obj=obj)

HSCimgLoader will automatically try to find a matching HSC object close to the coordinate. One can check whether it succeed by:

	>>> L.status

One can also change the image size in the same way as ``SDSSimgLoader``. 


Download stamps and psf
-----------------------
There are several options

	>>> status = L.make_stamps()
	>>> status = L.make_stamp(band='r')

	>>> status = L.make_psfs()
	>>> status = L.make_psf(band='r')

The ``status`` (bool) argument tells you whether the downloading was successful. 


Data Release
------------

The default data release is set to dr1 s16a_wide, one can specify the data release by:

	>>> L = HSCimgLoader(obj=obj, release_version='dr1', rerun='s15b_udeep')


Side note
---------

Instead of providing ``dir_obj`` to initiate the object, one can also provide ``dir_parent``. And directory under dir_parent with a name such as SDSSJ1234+4312 will be created under ``dir_parent`` as ``dir_obj``. It works for both obsobj and imgloaders. 

	>>> from bubbleimg.imgdownload import SDSSimgLoader
	>>>
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_parent = './'
	>>> L = SDSSimgLoader(ra=ra , dec=dec, dir_parent=dir_parent)
