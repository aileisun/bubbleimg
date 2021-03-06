*****
Batch
*****

``batch`` takes care of operations done to a sample of objects. It has three basic functions: ``build()``, ``iterlist()``, and ``compile_table()``. The batch class is still under construction. Please use ``hscBatch`` for now. The computationally intensive parts (``build()``, ``iterlist()``) are run in parallel using multiprocessing by default. 


Batch
=====

Create Batch
------------

To create a batch one needs to provide a directory path ``dir_batch`` (string), a sample catalog ``catalog`` (astropy table), and the survey of the sample ``survey`` (string). 

	>>> from bubbleimg.batch import Batch
	>>> from astropy.io import ascii
	>>> dir_batch = 'here/my_batch/'
	>>> lines = ['    ra       dec    ', 
				'--------- ----------', 
				'29.158592 -4.0001336', 
				'29.748938 -6.4643413', 
				'30.410525 -6.3772438', ]

	>>> catalog = ascii.read(lines)
	>>> survey = 'hsc'
	>>> b = Batch(dir_batch=dir_batch, catalog=catalog, survey=survey)

Alternatively, one can provide a path to a fits file that contains the catalog table instead of the catalog itself. 

	>>> fn_cat = 'catalog.fits'
	>>> b = Batch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey)


The catalog table should contain at least two columns ``ra``, ``dec``, representing the coordinates of the objects in degree (decimal J2000). 

The name of the batch will be assigned to be the last directory of ``dir_batch``. 

	>>> b.name
	'my_batch'

One could also provide ``dir_parent`` and ``name`` instead of ``dir_batch``. 

	>>> dir_parent = 'here/'
	>>> name = 'my_batch'
	>>> b = Batch(dir_parent=dir_parent, catalog=catalog, survey=survey, name=name)
	>>> b.name
	'my_batch'
	>>> b.dir_batch
	'here/my_batch/'


One can check the list of objects of batch. 

	>>> b.list
	<Table length=3>
	    ra       dec        obj_name   
	 float64   float64       str64     
	--------- ---------- --------------
	29.158592 -4.0001336 SDSSJ0156-0400
	29.748938 -6.4643413 SDSSJ0158-0627
	30.410525 -6.3772438 SDSSJ0201-0622


obj_naming_sys
--------------

By default the naming of the objects is SDSSJ plus the hhmmsddmm of the coordinates (sdss). The naming system can be changed using the argument ``obj_naming_sys``. 

	>>> b = Batch(dir_batch=dir_batch, fn_cat=fn_cat, survey=survey, obj_naming_sys='sdss')

Currently supported ``obj_naming_sys``: 
		- 'sdss'			: 'SDSSJ1000+1242'
		- 'sdss_precise'	: 'SDSSJ100013+124226'
		- 'j'				: 'J1000+1242'
		- 'j_precise'		: 'J100013+124226'

If the sample is too large such that there might be duplicated object names, consider using more precise naming systems. 


args_to_list
------------

One can have additional columns from the catalog to be included in list.csv by using ``args_to_list``.

	>>> lines = [
				'    ra       dec        z  ', 
				'--------- ---------- ------', 
				'29.748938 -6.4643413 0.4178', 
				'32.258361 -6.410809  0.4248', 
				'33.480903 -5.8281496 0.4434', 
				]

	>>> catalog = ascii.read(lines)
	>>> b = Batch(dir_batch=dir_batch, catalog=catalog, survey=survey, args_to_list=['z'])
	>>> b.list
	<Table length=3>
	    ra       dec        obj_name       z   
	 float64   float64       str64      float64
	--------- ---------- -------------- -------
	29.748938 -6.4643413 SDSSJ0158-0627  0.4178
	32.258361  -6.410809 SDSSJ0209-0624  0.4248
	33.480903 -5.8281496 SDSSJ0213-0549  0.4434



hscBatch
========

Create hscBatch
---------------

An ``hscBatch`` instance is created in the same way as ``Batch``. One does not need to set the ``survey`` argument and it's automatically set to ``'hsc'``. 

	>>> from bubbleimg.batch.hscbatch import hscBatch
	>>> dir_batch = 'here/my_batch/'
	>>> fn_cat = 'catalog.fits'
	>>> b = hscBatch(dir_batch=dir_batch, fn_cat=fn_cat)


Build hscBatch
--------------

Building the batch means to download the relevent data for each of the objects in the batch for the following analysis, and classify the objects into two classes: ``good`` and ``except``, where ``good`` objects have all the required files ready and ``except`` objects had problem getting all the files. 

One can do 

	>>> status = b.build()

By default, the ``bulid()`` function of hscBatch does the following steps to each of the objects: 

	- create an instance of hscimgLoader
		which tries to find an matching hsc counterpart 

	if that is successful then:
		- downlaod the 5 band stamp cutout images
		- download the 5 band psf
		- make a false color image of the stamp cutout
		- tries to find an sdss counterpart
		- download sdss spectrum


If the building was successful then the returning status is True. 
	>>> status
	True

If all the relevent files are successfully downloaded for an object, it will be inlcuded in the ``list_good`` list. (Here is just an example)

	>>> b.list_good
	<Table length=1>
	    ra       dec        obj_name   
	 float64   float64       str64     
	--------- ---------- --------------
	29.158592 -4.0001336 SDSSJ0156-0400

The object directories will be stored under dir_obj/good/. 

Otherwise it will be in the ``list_except`` list. 

	>>> b.list_except
	<Table length=2>
	    ra       dec        obj_name   
	 float64   float64       str64     
	--------- ---------- --------------
	29.748938 -6.4643413 SDSSJ0158-0627
	30.410525 -6.3772438 SDSSJ0201-0622

The object directories will be stored under dir_obj/except/. 


If you want to do the downloading again and overwrite the previously downloaded files. Do

	>>> status = b.build(overwrite=True)

If you want to include other columns from your input catalog into the batch list, include them as ``
listargs
--------


Customizing your build
----------------------

The default building setting is specified by the in ``self._func_build()`` of ``hscBatch``.

	>>> def _func_build(self, obj, overwrite=False, **kwargs):
	>>> 	"""
	>>> 	Params
	>>> 	------
	>>> 	obj
	>>> 	overwrite=False
	>>> 
	>>> 	**kwargs:
	>>> 		environment='iaa'
	>>> 
	>>> 	Return
	>>> 	------
	>>> 	status
	>>> 	"""
	>>> 
	>>> 	# setting
	>>> 	environment = kwargs.pop('environment', 'iaa')
	>>> 	humvi_bands = 'riz'
	>>> 
	>>> 	# running
	>>> 	L = imgdownload.hscimgLoader(obj=obj, environment=environment, **kwargs)
	>>> 
	>>> 	if L.status:
	>>> 		statuss = 	[ 
	>>> 					L.make_stamps(overwrite=overwrite), 
	>>> 					L.make_psfs(overwrite=overwrite), 
	>>> 					L.plot_colorimg(bands=humvi_bands, img_type='stamp', overwrite=overwrite)
	>>> 					L.add_obj_sdss(), 
	>>> 					L.obj.sdss.make_spec(overwrite=overwrite),
	>>> 					]
	>>> 
	>>> 		return all(statuss)
	>>> 	else:
	>>> 		return False


One can change how it's built by writing one's own ``func_build()``. This function has to take ``obj`` (instance of obsObj), see documentation for obsobj, and ``overwrite`` (bool), which specify whether to overwrite the downloaded files, as input arguments, and optionally other arguments as ``**kwargs``. This function has to return ``status`` to indicate whether the building of an object was successful. 

For example, one can define a very simple ``func_build()``,

	>>> def func_build(obj, overwrite=False):
	>>> 	"""
	>>> 	Params
	>>> 	------
	>>> 	obj
	>>> 	overwrite=False
	>>> 
	>>> 	Return
	>>> 	------
	>>> 	status
	>>> 	"""
	>>> 	status = obj.add_hsc()
	>>> 	
	>>> 	return status

which only tries to find an hsc counterpart and stores its info as hsc_xid.csv. 

To run it, one can do 
	>>> status = b.build(func_build)


iterlist
--------

Once a batch is built then you can perform operations, usually calculating values, on the built batch with ``iterlist()``. 

	>>> status = b.iterlist(func_iterlist, **kwargs)

By default the operation is applied to each of the ``good`` objects and the ``except`` objects will be ignored. 


You will need to define a function, for example, ``func_iterlist`` to be applied to each of the objects in the batch. For example:


	>>> def func_iterlist(obj, overwrite=False, **kwargs):
	>>> 	fn_testing = 'testing.txt'
	>>> 	fn = obj.dir_parent+fn_testing
	>>> 	print fn
	>>> 	with open(fn, 'a') as f:
	>>> 		f.write(obj.name+'\n')
	>>> 
	>>> 	return True
	>>>
	>>> status = b.iterlist(func_iterlist, **kwargs)
	>>> 

This function has to take ``obj`` (``obsobj`` instance) and ``overwrite`` (bool) as arguments, and optionally other arguments as ``**kwargs``. It should also return status (bool). 

Another example that queries the hsc forced catalog:


	>>> def func_iterlist_get_hsc_photoboj(obj, overwrite=False):
	>>> 	statuss = [
	>>> 				obj.add_hsc()
	>>> 				obj.hsc.load_photoobj(overwrite=overwrite)
	>>> 				]
	>>> 
	>>> 	return all(statuss)
	>>>
	>>> status = b.iterlist(func_iterlist, **kwargs)
	>>> 


This ``func_iterlist`` can take in additional arguments from the batch list, for example, to have redshift ``z`` as an additional argument one can have

	>>> def func_iterlist_make_spec_mag(obj, z, overwrite=False):
	>>> 	""" 
	>>> 	make file spec_mag.csv 
	>>> 
	>>> 	Params
	>>> 	------
	>>> 	obj
	>>> 	overwrite=False
	>>> 
	>>> 	Return
	>>> 	------
	>>> 	status
	>>> 	"""
	>>> 	s = bubbleimg.spector.Spector(obj=obj, z=z)
	>>> 
	>>> 	statuss = [
	>>> 				s.plot_spec(wfilters=True, wconti=True, overwrite=overwrite),
	>>> 				s.make_spec_mag(overwrite=overwrite),
	>>> 				]
	>>> 
	>>> 	return all(statuss)

Can call it by

	>>> statuss = b.iterlist(func_iterlist_make_spec_mag, listargs=['z'], overwrite=False)

This will create a one row table ``spec_mag.csv`` for each of the objects containing the spectroscopic magnitudes, which can be compiled over the entire sample by ``compile_table()``. The return value of ``iterlist()`` is a list containing the ``iterfunc()`` return value of each of the objects. 


If in a rare occasion where you want to iterate the function through the ``except`` list, do

	>>> status = b.iterlist(func_iterlist, listname='except', **kwargs)




compile_table
-------------

To compile the results table that resides in each of the object directories, one can do, for example

	>>> status = b.compile_table('spec_mag.csv', overwrite=True)

This will create ``dir_batch/spec_mag.csv`` that contains the ``spec_mag.csv`` for all of the objects in the list, including the ``exclude`` objects. Just that the content of the ``exclude`` object will be empty. 



steal_columns
-------------

If after building a batch you realized that there are columns in the catalog that you wanted to include as 'args_to_list' in list.csv (and all the compiled tables), you can make it up by using ``steal_columns()``. 

	>>> import astropy.table as at
	>>> cat = at.Table.read('catalog.fits')
	>>> b.steal_columns(tab=cat, colnames=['z'])

The list.csv, list_good.csv, and list_except.csv files will be overwrite and the new columns (in this example 'z') will appear. 

To have it in the compiled tables you will need to run ``compile_table`` again. 
	>>> status = b.compile_table('spec_mag.csv', overwrite=True)

If you want to undo ``steal_columns()``, use ``remove_columns()``

	>>> b.remove_columns(colnames=['z'])



Multiprocessin
==============

By the fault the computations in ``build()`` and ``iterlist()`` is run in parallel using multiprocessing, see <https://docs.python.org/2/library/multiprocessing.html>. More specifically, we uses the ``map`` function of ``multiprocessing.Pool``. The default number of processes is determined by the default of the map function. To change the number of processes, do: 

	>>> status = b.build(func_build, processes=2)

	or

	>>> statuss = b.iterlist(func_iterlist, processes=2)

To not use multiprocessing at all and run through the list using python for loop (!!!), set ``processes`` to -1. 

	>>> status = b.build(func_build, processes=-1)

	or

	>>> statuss = b.iterlist(func_iterlist, processes=-1)

