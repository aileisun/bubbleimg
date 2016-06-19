# do_batch.py
# ALS 2016/05/03

"""
once a batch is made (with bathc/list.txt and subdirectories), run certain
operations on the whole batch. objobs is not required.
"""
import os
from astropy.table import Table, hstack, vstack
import astropy.units as u

import smallfunc
reload(smallfunc)

import measureimg
reload(measureimg)

import denoiseimg
reload(denoiseimg)

import fitpsf


def do_batch(dir_batch, bandline=None, update=True):
    """
    Do a bunch of operations on a batch directory:

    0. make normalized OIII line maps
    1. calculate its noise level and compile batch noise table
    2. make denoised map
    3. make ISO measurements on the denoised map and compile the batch table
    4. join the measurement table with Mullaney table
    """
    isocut_rest = 3.e-15*u.Unit('erg / (arcsec2 cm2 s)')

    # #==== renoamlize
    # kwargs={'filename':'stamp-lOIII5008_I.fits','norm':1.e-15,'update':update}
    # do_mapjob_onbatch(dir_batch, smallfunc.dir_RenormalizeImg_fits,**kwargs)

    # kwargs={'filename':'stamp-conti-onOIIIscale_I.fits','norm':1.e-15,'update':update}
    # do_mapjob_onbatch(dir_batch, smallfunc.dir_RenormalizeImg_fits,**kwargs)

    # #==== calculate noise
    # kwargs={'filename':'stamp-lOIII5008_I_norm.fits','update':update}
    # do_mapjob_onbatch(dir_batch, denoiseimg.noiselevel.load_noiselevel,**kwargs)
    # # compile batch noise table
    # do_compiletable_onbatch(dir_batch,'noiselevel.csv')

    # # ==== do psf fitting
    # kwargs = dict(band=bandline, fileimg='stamp-lOIII5008_I_norm.fits',
    #               searchradius=5., fixb=True)
    # do_mapjob_onbatch(dir_batch, fitpsf.dir_fit_psf, **kwargs)


    # # #==== denoise
    # # kwargs={'filename':'stamp-lOIII5008_I_norm.fits','update':update}
    # # do_mapjob_onbatch(dir_batch, denoiseimg.dir_makedenoised_fits,**kwargs)
    # kwargs={'filename':'stamp-lOIII5008_I_norm_psfresidual.fits','update':update}
    # do_mapjob_onbatch(dir_batch, denoiseimg.dir_makedenoised_fits,**kwargs)


    # ==== deleting redundant files
    filenames=[
    "measureimg_iso_contoursdict_linblobmasked.pdf",
    "contoursdict_linblobmasked.pkl",
    "contoursdict_stamp-lOIII5008_I_norm_psfmodel_iso3e-15_alim0.pkl",
    "measureimg_iso_stamp-lOIII5008_I_norm_denoised_contours.pkl", 
    "measureimg_iso_stamp-lOIII5008_I_norm_denoised.ecsv", 
    "measureimg_iso_stamp-lOIII5008_I_norm_denoised.csv", 
    "measureimg_stamp-lOIII5008_I_norm_denoised.ecsv", 
    "measureimg_stamp-lOIII5008_I_norm_denoised.csv", ]
    kwargs=dict(filenames=filenames)
    do_mapjob_onbatch(dir_batch, smallfunc.dir_delete_files, **kwargs)


    # ==== make new iso measurements
    kwargs = {'isocut_rest': isocut_rest,
              'isoareallimit': 10, 'contrastr': 0.1,
              'update': update, 'toplot': True}
    do_mapjob_onbatch(dir_batch, measureimg.dir_doIsos, **kwargs)

    # # ==== join the measurement table with Mullaney table
    # smallfunc.joinmullaney(dir_batch,filename='measureimg_'+filename+'.ecsv')



def do_mapjob_onbatch(dir_batch, function, **kwargs):
    """
    apply function to all objects in the batch. The funciton takes parameter 
    dir_obj and optionally **kwargs, and give no returns. 

    data need to be organized in this way:
    dir_batch/
        list.txt <- contains a column of object 'OBJNAME'
        SDSSJ0012-0947/
        SDSSJ0033+0045/
        ...

    Parameters
    ------
    dir_batch: string
        directory path to batch
    function: function
        a function that takes parameter dir_obj

    Returns
    -----
    None

    Write Output
    -----
    depending on the function
    """
    print "doing map "+function.func_name

    listin=Table.read(dir_batch+'list.txt',format='ascii')
    listexclude=Table.read(dir_batch+'list_exclude.txt',format='ascii')

    for objname in listin['OBJNAME']:
        if objname not in listexclude['OBJNAME']:
            print objname
            dir_obj=dir_batch+objname+'/'
            function(dir_obj,**kwargs)


def do_reducejob_onbatch(dir_batch,function,filename_tabout,**kwargs):
    """
    apply function to all objects in the batch. The funciton takes parameter 
    dir_obj and optionally **kwargs, the function returns a table, and the 
    tables are compiled (each row marked by OBJNAME) and saved under 
    filename_tabout.

    Data need to be structured in this way:
    dir_batch/
        list.txt <- contains a column of object 'OBJNAME'
        SDSSJ0012-0947/
        SDSSJ0033+0045/
        ...

    Parameters
    ------
    dir_batch: string
        directory path to batch
    function: function
        a function that takes parameter dir_obj
        and return table, tableformat
    """
    # print "doing reduce "+function.func_name

    listin=Table.read(dir_batch+'list.txt',format='ascii')
    listexclude=Table.read(dir_batch+'list_exclude.txt',format='ascii')

    tabout=Table()
    tableformat='ascii.csv'

    for objname in listin['OBJNAME']:
        if objname not in listexclude['OBJNAME']:

            print objname
            dir_obj=dir_batch+objname+'/'

            # run function
            tabfunc, tableformat=function(dir_obj,**kwargs)

            # compile table
            tabobjname=Table([[objname]],names=['OBJNAME'])
            tabrow=hstack([tabobjname,tabfunc])
            tabout=vstack([tabout,tabrow])

    tabout.write(dir_batch+filename_tabout,format=tableformat)


def do_compiletable_onbatch(dir_batch,filename):
    """
    read tables in indivisual dir_obj under dir_batch and stack them to a 
    big table under dir_batch

    """
    def readtable(dir_obj,filename):
        # read format
        extension=os.path.splitext(filename)[1]
        if extension=='.csv': tableformat='ascii.csv'
        elif extension=='.ecsv': tableformat='ascii.ecsv'
        elif extension=='.fits': tableformat='fits'
        else: raise NameError('file extension of table unrecognized')

        return Table.read(dir_obj+filename,format=tableformat), tableformat

    print "compiling table "+filename
    kwargs={'filename':filename}
    do_reducejob_onbatch(dir_batch,readtable,filename,**kwargs)




# def main():
#     """
#     can be flexibly edited to perform desired functions 

#     now 
#     make denoised images for all objects in the batch
#     takes stamp-lOIII5008_I.fits and makes stamp-lOIII5008_I_denoised.fits
#     """

#     import denoiseimg
#     reload(denoiseimg)



#     dir_mullaney=external_links.dir_data_mullaney

#     batches=['mullaney_gi_T2','mullaney_gr_T2','mullaney_rz_T2','mullaney_ri_T2']
#     # batches=['mullaney_rz_T2']


#     for batch in batches:
#         print batch
#         dir_batch=dir_mullaney+batch+'/'

#         # # renoamlize
#         # kwargs={'norm':1.e-15}
#         # do_mapjob_onbatch(dir_batch, smallfunc.dir_RenormalizeImg_fits,**kwargs)

#         # # denoise
#         # kwargs={'filename':'stamp-lOIII5008_I_norm.fits'}
#         # do_mapjob_onbatch(dir_batch, denoiseimg.dir_makedenoised_fits,**kwargs)

#         # # recalculate noise
#         # kwargs={'filename':'stamp-lOIII5008_I_norm.fits','update':update}
#         # do_mapjob_onbatch(dir_batch, denoiseimg.noiselevel.load_noiselevel,**kwargs)

#         # # compile batch noise table
#         # do_compiletable_onbatch(dir_batch,'noiselevel.csv')

#         # filename='stamp-lOIII5008_I_norm'
#         filename='stamp-lOIII5008_I_norm_denoised'
#         # make iso measurements
#         kwargs={'filename':filename+'.fits',
#                 'isocut_rest':3.e-15*u.Unit('erg / (arcsec2 cm2 s)')}
#         do_mapjob_onbatch(dir_batch, measureimg.dir_MeasureImgIso,**kwargs)

#         # compile batch noise table
#         do_compiletable_onbatch(dir_batch,'measureimg_'+filename+'.csv')

        
#         smallfunc.joinmullaney(dir_batch,filename='measureimg_'+filename+'.csv')



# old code: 
    # #==== make iso measurements
    # filename='stamp-lOIII5008_I_norm_denoised'
    # kwargs={'fileimg':filename+'.fits',
    #         'isocut_rest': isocut_rest,
    #         'isoareallimit':10,
    #         'update':update}
    # do_mapjob_onbatch(dir_batch, measureimg.dir_MeasureImgIso, **kwargs)
    # # compile measurement table
    # # do_compiletable_onbatch(dir_batch,'measureimg_iso_'+filename+'.ecsv')
    # # do_compiletable_onbatch(dir_batch,'measureimg_iso_'+filename+'.csv')
