# smallfunc.py
# ALS 2016/05/04

"""
a temporary repository of small functions 
"""
import os
from astropy.io import fits
import astropy.table as at
import astropy.units as u

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

def joinmullaney(dir_batch, filein, fileout=None, filemullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/ALPAKA_extended.fits'):
    """
    join the table with mullaney table
    """
    import catalogue

    # write out
    if fileout is None:
        fileout = os.path.splitext(filein)[0]+'_joinmullaney'
 
    # files to join
    tabin=at.Table.read(dir_batch+filein,format='ascii.ecsv')
    tabmullaney=at.Table.read(filemullaney,format='fits')
    tabmullaney.rename_column('SDSSNAME','OBJNAME')

    tabmullaney_matched = catalogue.catalogue_util.selectRADEC(tabmullaney, tabin, radius=5., verbos=True)

    tabout=at.join(tabin, tabmullaney_matched, keys=['OBJNAME'], join_type='left')

    tabout.write(dir_batch+fileout+'.ecsv',format='ascii.ecsv')
    tabout.write(dir_batch+fileout+'.csv',format='ascii.csv')


def matchedmullaney(dir_batch, filein, fileout=None, filemullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/ALPAKA_extended.fits'):
    import catalogue
    reload(catalogue)
    # write out
    if fileout is None:
        fileout = os.path.splitext(filein)[0]+'_matchedmullaney'

    # files to join
    tabin = at.Table.read(dir_batch+filein, format='ascii.ecsv')
    tabmullaney = at.Table.read(filemullaney, format='fits')
    tabmullaney.rename_column('SDSSNAME', 'OBJNAME')

    tabout = catalogue.catalogue_util.selectRADEC(tabmullaney, tabin, radius=5., verbos=True)

    if len(tabout) != len(tabin):
        raise ValueError("len of tabout does not match with tabin")

    tabout.write(dir_batch+fileout+'.ecsv',format='ascii.ecsv')
    tabout.write(dir_batch+fileout+'.csv',format='ascii.csv')


def write_list_ran(dir_batch, filein):
    """
    write those three columns of the table as a new table
    """
    print "writing list_ran.csv"
    fileout = 'list_ran.csv'
    tabin = at.Table.read(dir_batch+filein, format='ascii')

    tabout = tabin['OBJNAME', 'RA', 'DEC']
        
    tabout.write(dir_batch+fileout, format='ascii.csv')

