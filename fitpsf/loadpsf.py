# loadpsf.py

"""
return sdss psf given dir_obj and band

download 'psField*.fit' file and translate that to psf file using readAtlasImages-v5_4_11 read_PSF as needed. 

the resulting psf are stored in, e.g. 'psf-r.fit', files under dir_obj
"""

import os
import astropy.table as at
from astropy.io import fits

def dir_load_psf(dir_obj, band):
    """
    load the psf img np.array for band 'band'. download/make files as nessisary. psField.fit file is deleted afterwards as it's too big. 
    """

    infile=dir_obj+'psf-'+band+'.fit'
    parentfile=dir_obj+'psField.fit'

    if not os.path.isfile(infile): # make needed files when they are absent
        if os.path.isfile(parentfile): 
            # parentfile exist, translate from parent
            dir_psField_to_psf(dir_obj, band=band)
        else:
            # parentfile doesn't exist, download from web, translate all bands
            dir_download_psField(dir_obj)
            for b in ['u','g','r','i','z']:
                dir_psField_to_psf(dir_obj, band=b)
            os.remove(parentfile)

    psf=fits.getdata(infile)
    return psf


def dir_psField_to_psf(dir_obj, band='r'):
    """ for a dir_obj translate psField to psf-(band).fit file given band    """
    codepath='/Users/aisun/Documents/Astro/codelab/readAtlasImages-v5_4_11/'
    infile=dir_obj+'psField.fit'
    outfile=dir_obj+'psf-'+band+'.fit'

    if os.path.isfile(infile):
        iband={'u':1,'g':2,'r':3,'i':4,'z':5,}
        ib=iband[band]

        # get row, col (position on image) from photoobj
        photoobj = at.Table.read(dir_obj+'PhotoObj.csv',format='ascii.csv')
        row=photoobj['rowc_'+band]
        col=photoobj['colc_'+band]

        # running command
        command=codepath+'read_PSF '+infile+' '+'%d'%ib+' '+'%.4f'%row+' '+'%.4f'%col+' '+outfile
        os.system(command)
    else: 
        raise ValueError('psField.fit does not exist')


def dir_download_psField(dir_obj):
    """ for a dir_obj read in xid and download psField"""
    xid = at.Table.read(dir_obj+'xid.csv',format='ascii.csv')
    download_psField(xid, dir_out=dir_obj, filenameout='psField.fit')


def download_psField(xid, dir_out='./', filenameout=None):
    """
    download the psField to directory dir_out as 'psField.fit' unless otherwise specified. 
    """
    import urllib
    print "downloading psField file for psf"

    def get_urlpsField(run, rerun, camcol, field):
        urlsdssbase='http://data.sdss3.org/sas/dr12/'
        urlpsField='boss/photo/redux/'+'%d'%rerun+'/'+'%d'%run+'/objcs/'+'%d'%camcol+'/psField-'+'%06d'%run+'-'+'%d'%camcol+'-'+'%04d'%field+'.fit'
        return urlsdssbase+urlpsField

    # sanity check
    if not isinstance(xid, at.Table) or len(xid)!=1:
        raise ValueError("xid is not a one row table as needed")

    # constructing URL
    run, rerun, camcol, field=[xid[key][0] for key in ['run', 'rerun', 'camcol', 'field']]
    url=get_urlpsField(run, rerun, camcol, field)

    # output setup
    if filenameout is not None: filenameout='psField.fit'
    else: filenameout=url.rsplit('/', 1)[-1]

    # download
    testfile = urllib.URLopener()
    testfile.retrieve(url, dir_out+filenameout)

