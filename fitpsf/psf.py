# psf.py

"""
class psf 
"""


class psf(object):

    def __init__(self, dir_obj, band):
        """
        Params
        ----------
        dir_obj: string
            directory of the object, e.g., 'J1356+1026'
        band: string
            band of the psf, e.g., 'r'
        """
        self.dir_obj = dir_obj
        self.band = band

    def get_psfimg(self):
        filename = self.dir_obj+'psf-'+self.band+'.fit'

        if not os.path.isfile(filename): 
            self.make_psf_fits()

        psfimg = fits.getdata(filename)
        return psfimg

    def make_psf_fits(self):

        fileout = self.dir_obj+'psf-'+self.band+'.fit'
        filein = self.dir_obj+'psField.fit'
        ....... to continue ..........  

        # make needed files when they are absent
            if os.path.isfile(parentfile): 
                # parentfile exist, translate from parent
                dir_psField_to_psf(dir_obj, band=band)
            else:
                # parentfile doesn't exist, download from web, translate all bands
                dir_download_psField(dir_obj)
                for b in ['u','g','r','i','z']:
                    dir_psField_to_psf(dir_obj, band=b)
                os.remove(parentfile)


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

