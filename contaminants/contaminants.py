# contaminants.py
# ALS 2016/06/20

"""
for each object find out number of potential contaminants 
"""
import os
import numpy as np
import astropy.table as at
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
import astropy.units as u

def dir_reduce_contaminantcount(dir_obj, bandmag='r', min_modelMags=[18, 21], max_seps=[3, 6, 18]):
    """
    summarize the contaminants table for do_batch.do_reducejob_onbatch to compile a big table
    """
    print "running dir_reduce_contaminantcount"
    filein = dir_obj+'contaminants.csv'
    tabin = at.Table.read(filein, format='ascii.csv')

    min_modelMags = np.sort(min_modelMags)
    tabout = at.Table()

    for i, min_modelMag in enumerate(min_modelMags):
        tag_mag = '_mag'+bandmag+str(min_modelMag)
        sel_mag = tabin['modelMag_'+bandmag] < min_modelMag

        for j, max_sep in enumerate(max_seps): 
            tag_sep = '_sep'+str(max_sep)
            sel_sep = tabin['seperation'] < max_sep

            # select
            tabsel = tabin[sel_mag & sel_sep]
            # count n
            ncont = len(tabsel)
            # assemble types
            conttypes_list = [str(row['type']) for row in tabsel]
            conttypes_string = ''.join(np.sort(conttypes_list))

            tabnew = at.Table([[ncont], [conttypes_string]], names=['n_contaminants'+tag_mag+tag_sep, 'types'+tag_mag+tag_sep])

            tabout = at.hstack([tabout, tabnew])

    return tabout, 'ascii.csv'



def dir_find_contaminants(dir_obj, bandmag='r', min_modelMag=21, radius=18*u.arcsec, update=True):
    """
    for dir_obj find contaminatns and write contaminants.csv

    Params
    ------
    dir_obj
    bandmag='r'
        which band to do mag selection on
    min_modelMag=21
        the faintest obj to be considered as contaminants
    radius=18*u.arcsec
        search radius

    Write Output
    ------
    contaminants.csv: table with columns:
          ra            dec       mode  type modelMag_u modelMag_g modelMag_r modelMag_i modelMag_z   seperation 
    """
    # define filenames
    filexid = dir_obj+'xid.csv'
    fileout = dir_obj+'contaminants.csv'

    if not os.path.isfile(fileout) or update:
        # get xid
        xid = at.Table.read(filexid, format='ascii.csv', comment='#')
        ra, dec = xid['ra', 'dec'][0]

        # find contaminants
        kwargs = dict(bandmag=bandmag, min_modelMag=min_modelMag, radius=radius)
        tabcon = find_contaminants(ra, dec, **kwargs)

        tabcon.write(fileout, format='ascii.csv')
    else: 
        print "skip dir_find_contaminants as file exists"


def find_contaminants(ra, dec, bandmag='r', min_modelMag=21, radius=25*u.arcsec):
    """
    Find contaminants around coordinate. 

    Returns
    ------
    results: at table
        with cols: 
      ra            dec       mode  type modelMag_u modelMag_g modelMag_r modelMag_i modelMag_z   seperation 


    """
    photoobj_fields=['ra', 'dec', 'mode', 'type', 'modelMag_u',  'modelMag_g',  'modelMag_r',  'modelMag_i',  'modelMag_z',]

    c = SkyCoord(ra, dec, 'icrs', unit='deg')

    # find the galaxy itself
    result_self = SDSS.query_region(c, radius=3.*u.arcsec, spectro=False, photoobj_fields=photoobj_fields)
    result_self = result_self[result_self['mode'] == 1]  # pick primay
    result_self = result_self[np.any([result_self['type'] == 3, result_self['type'] == 6], axis=0)]  # pick galaxy or star

    if len(result_self) == 0: 
        raise ValueError("cannot find object itself")
    elif len(result_self) > 1: 
        print "more than one obj found, assigning the closest one as self"
        # likely merging galaxies
        separation=np.zeros(len(result_self), dtype=bool)
        for i, row in enumerate(result_self):
            c2 = SkyCoord(row['ra'], row['dec'], 'icrs', unit='deg')
            separation[i] = c.separation(c2).to(u.arcsec).value
        result_self = result_self[np.argmin(separation)]
    if result_self['type'] == 6: 
        print "obj self is star"

    # find all the other stuffs
    result = SDSS.query_region(c, radius=radius, spectro=False, photoobj_fields=photoobj_fields)

    # select bright primary contaminants
    result = result[result['mode']==1]  # pick primay
    result = result[result['modelMag_'+bandmag]<min_modelMag]  # pick bright 

    # remove itself from contaminant
    iself = np.where(result == result_self)[0]
    result.remove_rows(iself)

    # add seperation
    if len(result) == 0:
        col = at.Column(name='seperation', dtype=float)
        result.add_column(col)
    else:
        result['seperation'] = 0.
        for i, row in enumerate(result):
            c2 = SkyCoord(row['ra'], row['dec'], 'icrs', unit='deg')
            result['seperation'][i] = c.separation(c2).to(u.arcsec).value

    return result
