# smallfunc.py
# ALS 2016/05/04

"""
a temporary repository of small functions 
"""
import os
import numpy as np

from astropy.io import fits
import astropy.table as at
import astropy.units as u

import external_links

def dir_RenormalizeImg_fits(dir_obj,filename='stamp-lOIII5008_I.fits',norm=1.e-15,update=False):
    """
    make 'stamp-lOIII5008_I_norm.fits' which is scaled up by 1.e15
    """

    filein=dir_obj+filename
    fileout=dir_obj+os.path.splitext(filename)[0]+'_norm.fits'

    if not os.path.isfile(fileout) or update:
        # read in 
        img=fits.getdata(filein)
        header=fits.getheader(filein)

        img_new=img/norm
        header['HISTORY']="Renormalized by a factor of 1 over "+str(norm)
        header['BUNIT']='%.1e'%norm+" "+(header['BUNIT'])

        prihdu = fits.PrimaryHDU(img_new, header=header)
        prihdu.writeto(fileout, clobber=True)
    else: 
        print "skipping dir_RenormalizeImg_fits as file exists"


def dir_delete_file(dir_obj, filename):
    """ delete file "filename" in dir_obj """
    if os.path.isfile(dir_obj+filename):
        print "deleting file "+filename
        os.remove(dir_obj+filename)
    else:
        print "skip deleting file "+filename


def dir_delete_files(dir_obj, filenames):
    """ delete file "filename" in dir_obj """
    for filename in filenames: 
        dir_delete_file(dir_obj, filename)


def dir_fix_units_measureimg_iso(dir_obj):
    """ to fix a bug """ 
    pixelsize=0.396
    pixunit=u.Unit('arcsec')
    pixscale= pixelsize * (pixunit)

    filename = "measureimg_iso"
    tab = at.Table.read(dir_obj+filename+".ecsv", format='ascii.ecsv')
    if tab['isoshape_dferetper'].unit != pixunit and tab['isoshape_theta_dferetper'].unit != u.Unit('degree'):
        tab['isoshape_dferetper'] = tab['isoshape_dferetper']*pixelsize
        tab['isoshape_dferetper'].unit = pixunit
        tab['isoshape_theta_dferetper'].unit = u.Unit('degree')

        tab.write(dir_obj+filename+".ecsv", format='ascii.ecsv')
        tab.write(dir_obj+filename+".csv", format='ascii.csv')
    else: 
        print "skip dir_fix_units_measureimg_iso as units are correct"


def batch_joinmullaney(dir_batch, filein, fileout=None, filemullaney=None):
    """
    join the table with mullaney table
    """
    import catalogue

    # write out
    if fileout is None:
        fileout = os.path.splitext(filein)[0]+'_joinmullaney'
    if filemullaney is None:
        filemullaney = external_links.file_mullaney_extended
 
    # files to join
    tabin=at.Table.read(dir_batch+filein,format='ascii')
    tabmullaney=at.Table.read(filemullaney,format='fits')
    tabmullaney.rename_column('SDSSNAME','OBJNAME')

    tabmullaney_matched = catalogue.catalogue_util.selectRADEC(tabmullaney, tabin, radius=5., verbos=True)

    tabout=at.join(tabin, tabmullaney_matched, keys=['OBJNAME'], join_type='left')

    tabout.write(dir_batch+fileout+'.ecsv',format='ascii.ecsv')
    tabout.write(dir_batch+fileout+'.csv',format='ascii.csv')


def batch_matchedmullaney(dir_batch, filein, fileout=None, filemullaney=None):
    import catalogue
    reload(catalogue)
    # write out
    if fileout is None:
        fileout = os.path.splitext(filein)[0]+'_matchedmullaney'
    if filemullaney is None:
        filemullaney = external_links.file_mullaney_extended

    # files to join
    tabin = at.Table.read(dir_batch+filein, format='ascii')
    tabmullaney = at.Table.read(filemullaney, format='fits')
    tabmullaney.rename_column('SDSSNAME', 'OBJNAME')

    tabout = catalogue.catalogue_util.selectRADEC(tabmullaney, tabin, radius=5., verbos=True)

    if len(tabout) != len(tabin):
        raise ValueError("len of tabout does not match with tabin")

    tabout.write(dir_batch+fileout+'.ecsv',format='ascii.ecsv')
    tabout.write(dir_batch+fileout+'.csv',format='ascii.csv')


def batch_matchWISE(dir_batch, filein, fileout=None):
    """
    for a given table add new WISE columns
    """
    import catalogue

    if fileout is None:
        fileout = filein.split('.')[0]+'_WISE'
    tabin = at.Table.read(dir_batch+filein, format='ascii.ecsv')

    catalogue.crossmatch_util.matchWISE(tabin)
    # catalogue.crossmatch_util.matchIRAS(tabin)

    # adding MIR, FIR infered quantities
    catalogue.conversion_util.inferALLinTable(tabin)

    tabin.write(dir_batch+fileout+'.csv', format='ascii.csv')
    tabin.write(dir_batch+fileout+'.ecsv', format='ascii.ecsv')


def batch_write_list_ran(dir_batch, filein, fileout = 'list_ran.csv'):
    """
    write those three columns of the table as a new table
    """
    print "writing list_ran.csv"
    tabin = at.Table.read(dir_batch+filein, format='ascii')

    tabout = tabin['OBJNAME', 'RA', 'DEC']
        
    tabout.write(dir_batch+fileout, format='ascii.csv')


def batch_compile_results_physical(dir_batch):
    """
    compile a big table of relevent results in physical units

    Params
    ------
    dir_batch: string
        where the files locate
    rferet_bias: float
        the bias of rferet in arcsec -> to create rferets_debiased

    """
    print "compiling results physical"
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # out
    file_out = 'results_phys'
    # 
    file_list = 'list_ran.csv'
    file_mull = 'measureimg_iso_matchedmullaney_WISE.ecsv'
    # 
    file_blob_psfr = 'measureimg_iso_contours_blob_psfresid.ecsv'
    file_blob_psfr_m = 'measureimg_iso_contours_blob_psfresid_galmasked_psfmasked.ecsv'
    file_blob_psfrc = 'measureimg_iso_contours_blob_psfresidctr.ecsv'
    file_blob_psfrc_m = 'measureimg_iso_contours_blob_psfresidctr_galmasked_psfmasked.ecsv'
    file_gmaskc = 'measureimg_iso_contours_galmaskctr.ecsv'


    # 
    tab_list = at.Table.read(dir_batch+file_list, format='ascii.csv')
    tab_mull = at.Table.read(dir_batch+file_mull, format='ascii.ecsv')
    # 
    tab_blob_psfr = at.Table.read(dir_batch+file_blob_psfr, format='ascii.ecsv')
    tab_blob_psfr_m = at.Table.read(dir_batch+file_blob_psfr_m, format='ascii.ecsv')
    tab_blob_psfrc = at.Table.read(dir_batch+file_blob_psfrc, format='ascii.ecsv')
    tab_blob_psfrc_m = at.Table.read(dir_batch+file_blob_psfrc_m, format='ascii.ecsv')
    tab_gmaskc = at.Table.read(dir_batch+file_gmaskc, format='ascii.ecsv')


    z = tab_mull['Z']
    kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(z)

    tab_z = at.Table([z, kpc_per_arcsec], names=['Z', 'kpc_per_arcsec'])

    col_list = ['OBJNAME', 'RA', 'DEC']
    col_mull = ['AGN_TYPE', 'SDSS_OIII_EW', 'SDSS_OIII_EWERR', 'OIII_5007_FWHM', 'OIII_5007_FWHM_ERR', 'OIII_5007B_FWHM', 'OIII_5007B_FWHM_ERR', 'OIII_5007T_LUM', 'OIII_5007_LUM_ERR', 'OIII_5007B_LUM_ERR', 'OIII_5007AVG_FWHM', 'Lbol_OIII', 'WISE_nuLnu_rf8um', 'WISE_nuLnu_rf15um', 'WISE_nuLnu_rf22um', 'Lbol_8um', 'Lbol_15um', 'Lbol_22um']

    # hstacking list and mullaney
    tab_out = at.hstack([tab_list[col_list], tab_z, tab_mull[col_mull]])

    # hstacking isos

    tabs = [tab_blob_psfr, tab_blob_psfr_m, tab_blob_psfrc, tab_blob_psfrc_m, tab_gmaskc]
    tag_imgs = ['blobpsfr', 'blobpsfrm', 'blobpsfrc', 'blobpsfrcm', 'gmaskc']
    for i in range(len(tabs)):
        # assigning
        tab = tabs[i]
        tag_img = tag_imgs[i]
        # conversion
        areas = tab['isoshape_area']*kpc_per_arcsec**2
        rferets = tab['isoshape_rferetmax']*kpc_per_arcsec
        rferets_debiased =debias_rferetmax(tab['isoshape_rferetmax'])*kpc_per_arcsec
        dferets = tab['isoshape_dferetmax']*kpc_per_arcsec
        aspectrs = tab['isoshape_aspectr']
        # assemble
        tab_isophys = at.Table(data=[areas, rferets, rferets_debiased, dferets, aspectrs], names=[tag_img+'_area', tag_img+'_rferetmax', tag_img+'_rferetmax_debiased', tag_img+'_dferetmax', tag_img+'_aspectr'])
        tab_out = at.hstack([tab_out, tab_isophys])

    tab_out.write(dir_batch+file_out+'.ecsv', format='ascii.ecsv')
    tab_out.write(dir_batch+file_out+'.csv', format='ascii.csv')


def debias_rferetmax(r_biased):
    """ 
    according to the bias of the denoised rferetmax measurement concluded from the simulation, we correct the rferetmax measurements. 

    r are in arcsecs
    """
    def correction(r_biased):
        if r_biased < 3.5: 
            r_unbiased = r_biased*0.571
        else: 
            r_unbiased = r_biased-1.5
        return r_unbiased

    if isinstance(r_biased, float):
        r_unbiased = correction(r_biased)
    else: 
        r_unbiased = np.array([correction(r) for r in r_biased])

    return r_unbiased