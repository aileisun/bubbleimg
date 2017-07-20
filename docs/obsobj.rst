******************
Observation Object
******************

Observation object ``obsobj`` is one basic element of bubbleimg. It represents an object (basically a galaxy), and contains functionalities to find matching counterparts in other data bases. 


Create obsobj
=============

It is defined by three paramters ``ra``, ``dec``, and ``dir_obj``, where the coordinates are in decimal degree (J2000), and ``dir_obj`` is a string of the path to a directory where the data for this object will be stored. 

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

Each object will be assigned a name according to its coordinate. 

	>>> obj.name
	'SDSSJ1000+1242'

SDSS
====

Find the closest counterpart
----------------------------

To find SDSS counterpart:

	>>> status = obj.add_sdss()

If status is True then a counterpart is found. If more than one object is within the search radius then only the closest one will be identified. Currently it only looks for spectroscopy objects. Photometric only objects are ignored. The relevent info will be stored as dir_obj+'sdss_xid.csv'. If successful one can access the info of the SDSS object. 

To get the entire xid table which contains run, rerun, camcol, field, as well as, mjd, fiberID, specobjid, instrument, etc: 

	>>> obj.sdss.xid
		<Table length=1>
	      ra          dec             objid         run  ... run2d instrument sciencePrimary
	   float64      float64           int64        int64 ... int64    str4        int64     
	------------- ------------ ------------------- ----- ... ----- ---------- --------------
	150.054773552 12.707302737 1237664106315579464  4338 ...    26       SDSS              1

One can also access only the desired valuable only by, for example:

	>>> obj.sdss.camcol
	4


Specify Searching Parameters
----------------------------

One can specify the data release (default dr13) by
	
	>>> status = obj.add_sdss(data_release='dr13')

And the search radius (default 2 arcsec) by

	>>> status = obj.add_sdss(search_radius=2.*u.arcsec)

One can load the entire photoobj table with the option

	>>> status = obj.add_sdss(toload_photoobj=True)
	>>> obj.sdss.photoobj
	<Table length=1>
	       objID        skyVersion  run  rerun ...    TAI_r         TAI_i         TAI_z    
	       int64          int64    int64 int64 ...   float64       float64       float64   
	------------------- ---------- ----- ----- ... ------------ ------------- -------------
	1237664106315579464          2  4338   301 ... 4578717964.0 4578718035.74 4578718179.16	


One can load the spectrum by:

	>>> sp = obj.sdss.get_spec()

obsobj will first try to look for the relevent files locally, e.g. dir_obj+'sdss_xid.csv'. If it doesn't exist then do the query online. 


HSC
===

STARs account
-------------

One needs to provide the STARs username and password to access HSC data. The easiest way is to set your STARS username and passwords as environmental variable. In the command line:

	>>> export HSC_SSP_CAS_USERNAME
	>>> read -s HSC_SSP_CAS_USERNAME
	>>> export HSC_SSP_CAS_PASSWORD
	>>> read -s HSC_SSP_CAS_PASSWORD

Otherwise you will be asked to enter it manually each time bubbleimg accesses HSC data. 


Find the closest counterpart
----------------------------

Similarly for HSC, one can do :

	>>> status = obj.add_hsc()

If status is True then a counterpart is found. Similarly, only the closest one within the search radius is saved. 

One can access the basic ids of the object by: 

	>>> obj.hsc.xid

or

	>>> obj.hsc.tract
	>>> obj.hsc.patch


Specify Searching Parameters
----------------------------

The default setting of the search is set to the following. 

	>>> import astropy.units as u
	>>> status = obj.add_hsc(rerun='s16a_wide', release_version='dr1', search_radius=2.*u.arcsec)

Please change the parameters accordingly. 

