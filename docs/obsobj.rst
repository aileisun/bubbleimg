******************
Observation Object
******************

Observation object ``obsobj`` is one basic element of bubbleimg. It represents an object (basically a galaxy), and contains functionalities to find matching counterparts in other data bases. 


Define obsobj
=============

It is defined by three paramters ``ra``, ``dec``, and ``dir_obj``, where the coordinates are in degree decimal J2000 system, and ``dir_obj`` is a string of the path to a directory where the data for this object will be stored. 

	>>> from bubbleimg import obsobj
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_obj = './here/'
	>>> obj = obsobj.obsObj(ra=ra, dec=dec, dir_obj=dir_obj)

Alternatively, one can provide ``dir_parent`` instead of ``dir_obj``. A directory under dir_parent with a name such as SDSSJ1234+4312 will be created under ``dir_parent`` as ``dir_obj``. 

	>>> from bubbleimg import obsobj
	>>> ra = 150.0547735
	>>> dec = 12.7073027
	>>> dir_parent = './'
	>>> obj = obsobj.obsObj(ra=ra, dec=dec, dir_parent=dir_parent)


Find SDSS Counterparts
======================

Currently only two data bases are supported -- SDSS and HSC. One can do

	>>> status = obj.add_sdss()

If status is True and a counterpart is found. Currently it only looks for spectroscopy objects. Photometric only objects are ignored. The relevent info will be stored as dir_obj+'sdss_xid.csv'. If successful one can access the info by, for example:

	>>> obj.xid
	>>> obj.camcal

One can specify the data release (default dr13) by
	
	>>> status = obj.add_sdss(data_release='dr13')

And the search radius (default 2 arcsec) by

	>>> status = obj.add_sdss(search_radius=2.*u.arcsec)

One can load the entire photoobj table with the option

	>>> status = obj.add_sdss(toload_photoobj=True)
	>>> obj.photobj


One can load the spectrum by:

	>>> sp = obj.get_spec()

obsobj will first try to look for the relevent files locally, e.g. dir_obj+'sdss_xid.csv'. If it doesn't exist then do the query online. 



Find HSC counterparts
=====================
...