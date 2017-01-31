# do_batch.py
# ALS 2016/05/03

"""
once a batch is made (with bathc/list.txt and subdirectories), run certain
operations on the whole batch. objobs is not required.
"""
import os
import astropy.table as at
import astropy.units as u

import smallfunc
reload(smallfunc)

import measureimg
reload(measureimg)

import denoiseimg
reload(denoiseimg)

import contaminants
reload(contaminants)

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

    #==== renoamlize
    kwargs={'filename':'stamp-lOIII5008_I.fits','norm':1.e-15,'update':update}
    do_mapjob_onbatch(dir_batch, smallfunc.dir_RenormalizeImg_fits,**kwargs)

    kwargs={'filename':'stamp-conti-onOIIIscale_I.fits','norm':1.e-15,'update':update}
    do_mapjob_onbatch(dir_batch, smallfunc.dir_RenormalizeImg_fits,**kwargs)

    #==== calculate noise
    kwargs={'filename':'stamp-lOIII5008_I_norm.fits','update':update}
    do_mapjob_onbatch(dir_batch, denoiseimg.noiselevel.load_noiselevel,**kwargs)
    # compile batch noise table
    do_compiletable_onbatch(dir_batch, filein='noiselevel.csv')

    # ==== do psf fitting
    kwargs = dict(band=bandline, fileimg='stamp-lOIII5008_I_norm.fits',
                  searchradius=5., fixb=True)
    do_mapjob_onbatch(dir_batch, fitpsf.dir_fit_psf, **kwargs)


    # #==== denoise
    kwargs={'filename':'stamp-lOIII5008_I_norm.fits','update':update}
    do_mapjob_onbatch(dir_batch, denoiseimg.dir_makedenoised_fits,**kwargs)
    kwargs={'filename':'stamp-lOIII5008_I_norm_psfresidual.fits','update':update}
    do_mapjob_onbatch(dir_batch, denoiseimg.dir_makedenoised_fits,**kwargs)

    delete_redundant_files(dir_batch)

    # ==== make new iso measurements
    kwargs = {'isocut_rest': isocut_rest,
              'isoareallimit': 10, 'contrastr': 0.1,
              'update': update, 'toplot': True, 'toclean': False}
    do_mapjob_onbatch(dir_batch, measureimg.dir_doIsos, **kwargs)


    # ==== compile tables
    kwargs = dict(filein='measureimg_iso.ecsv',
                  selectkey='filecontoursdict',
                  selectvalues=['contours_blob.pkl',
                                'contours_blob_psfresid.pkl',
                                'contours_blob_psfresid_galmasked.pkl',
                                'contours_blob_psfresid_galmasked_psfmasked.pkl',
                                'contours_blobctr.pkl',
                                'contours_blob_psfresidctr.pkl',
                                'contours_blob_psfresidctr_galmasked.pkl',
                                'contours_blob_psfresidctr_galmasked_psfmasked.pkl',
                                'contours_galctr.pkl',
                                'contours_psf.pkl',
                                'contours_galmask.pkl',
                                'contours_galmaskctr.pkl',
                                'contours_psfmask.pkl',],)
    do_compiletables_onbatch(dir_batch, **kwargs)

    # # # ==== join the measurement table with Mullaney table
    # kwargs = dict(filein='measureimg_iso_contours_blob.ecsv', 
    #               fileout='measureimg_iso_matchedmullaney')
    # smallfunc.batch_matchedmullaney(dir_batch, **kwargs)

    # # match with wise
    # kwargs = dict(filein='measureimg_iso_matchedmullaney.ecsv', 
    #               fileout='measureimg_iso_matchedmullaney_WISE')

    # smallfunc.batch_matchWISE(dir_batch, **kwargs)

    # ==== write list_final
    kwargs = dict(filein='measureimg_iso_contours_blob.ecsv')
    smallfunc.batch_write_list_ran(dir_batch, **kwargs)

    # ==== find contaminants
    kwargs = dict(bandmag='r', min_modelMag=21, radius=25*u.arcsec, update=False)
    do_mapjob_onbatch(dir_batch, contaminants.dir_find_contaminants, **kwargs)

    # compile big contaminants.csv table
    fileout = 'contaminantcount.csv'
    do_reducejob_onbatch(dir_batch, contaminants.dir_reduce_contaminantcount, fileout=fileout)

    # smallfunc.batch_compile_results_physical(dir_batch)

def delete_redundant_files(dir_batch):
    """ deleting redundant files from old versions """
    filenames=[
    "measureimg_iso_contoursdict_linblobmasked.pdf",
    "measureimg_iso_stamp-lOIII5008_I_norm_denoised.pdf",
    "stamp-lOIII5008_I_norm_denoised_iso.pdf",
    "contoursdict_linblobmasked.pkl",
    "contoursdict_stamp-lOIII5008_I_norm_psfmodel_iso3e-15_alim0.pkl",
    "measureimg_iso_stamp-lOIII5008_I_norm_denoised_contours.pkl", 
    "measureimg_iso_stamp-lOIII5008_I_norm_denoised.ecsv", 
    "measureimg_iso_stamp-lOIII5008_I_norm_denoised.csv", 
    "measureimg_stamp-lOIII5008_I_norm_denoised.ecsv", 
    "measureimg_stamp-lOIII5008_I_norm_denoised.csv", 
    "measureimg_iso.csv",
    "measureimg_iso.ecsv",]
    kwargs=dict(filenames=filenames)
    do_mapjob_onbatch(dir_batch, smallfunc.dir_delete_files, **kwargs)


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


    listin=at.Table.read(dir_batch+'list.txt',format='ascii')

    filepath_listexclude = dir_batch+'list_exclude.txt'
    if os.path.isfile(filepath_listexclude):
        listexclude = at.Table.read(filepath_listexclude,format='ascii')
    else: 
        print "Ignoring listexclude"
        listexclude = at.Table([[]], names=['OBJNAME'])

    for objname in listin['OBJNAME']:
        if objname not in listexclude['OBJNAME']:
            print objname
            dir_obj = dir_batch+objname+'/'
            function(dir_obj,**kwargs)


def do_compiletables_onbatch(dir_batch, filein, selectkey, selectvalues):
    """
    for each selectvalue compile a big table with its name

    Params
    -----
    filein: str
        filename of the obj table to compile
    selectkey: str
        within the obj table based on what col to split rows
        e.g., filecontoursdict
    selectvalues: list
        a list of values for the key

    """
    for selectvalue in selectvalues:
        inprefix = os.path.splitext(filein)[0]
        insuffix = os.path.splitext(filein)[1] # file extension
        outprefix = os.path.splitext(str(selectvalue))[0]
        fileout = inprefix+'_'+outprefix+insuffix
        do_compiletable_onbatch(dir_batch, filein, fileout, selectkey=selectkey, selectvalue=selectvalue)        


def do_compiletable_onbatch(dir_batch, filein, fileout='', selectkey=None, selectvalue=None):
    """
    read tables 'filein' in indivisual dir_obj under dir_batch and stack them 
    to a big table 'fileout' under dir_batch. 
    """
    if fileout == '': 
        fileout = filein
    print "compiling table "+fileout+" from "+filein
    kwargs = {'filename': filein, 'selectkey': selectkey, 'selectvalue': selectvalue}
    do_reducejob_onbatch(dir_batch, readtable, fileout=fileout, **kwargs)


def readtable(dir_obj, filename, selectkey=None, selectvalue=None):
    """ read table with format """
    # read extension
    extension = os.path.splitext(filename)[1]
    formatdict = {'.csv': 'ascii.csv', '.ecsv': 'ascii.ecsv', '.fits': 'fits'}
    tableformat = formatdict[extension]

    # read table
    tab = at.Table.read(dir_obj+filename, format=tableformat)

    # select table rows
    if selectkey is not None:
        tab = tab[tab[selectkey] == selectvalue]

    if len(tab) > 1:
        print "Warning: multiple entries per object being compiled"

    return tab, tableformat


def do_reducejob_onbatch(dir_batch, function, fileout, **kwargs):
    """
    apply function to all objects in the batch. The funciton takes parameter 
    dir_obj and optionally **kwargs, the function returns a table, and the 
    tables are compiled (each row marked by OBJNAME) and saved under 
    fileout.

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

    listin=at.Table.read(dir_batch+'list.txt',format='ascii')

    filepath_listexclude = dir_batch+'list_exclude.txt'
    if os.path.isfile(filepath_listexclude):
        listexclude = at.Table.read(filepath_listexclude,format='ascii')
    else: 
        print "Ignoring listexclude"
        listexclude = at.Table([[]], names=['OBJNAME'])

    tabout=at.Table()
    tableformat='ascii.csv'

    for i, objname in enumerate(listin['OBJNAME']):
        if objname not in listexclude['OBJNAME']:

            # printing
            print objname

            dir_obj = dir_batch+objname+'/'

            # run function
            tabfunc, tableformat = function(dir_obj, **kwargs)

            # compile table
            ra, dec = listin[i]['RA'], listin[i]['DEC']
            tabheader = at.Table([[objname], [ra], [dec]], names=['OBJNAME', 'RA', 'DEC'])
            tabrow = at.hstack([tabheader, tabfunc])
            tabout = at.vstack([tabout, tabrow])

    tabout.write(dir_batch+fileout, format=tableformat, overwrite=True)
    tabout.write(dir_batch+fileout.split('.')[0]+'.csv', format='ascii.csv', overwrite=True)


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


    # # ==== classify objects
    # filein = 'measureimg_iso_contours_blob.csv'
    # smallfunc.write_list_ran(dir_batch, filein)
    # classify.batch_classify(dir_batch)
