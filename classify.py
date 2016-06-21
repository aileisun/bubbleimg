# classify.py

import numpy as np
import astropy.table as at

def batch_classify(dir_batch):
    """ for a given batch, run classifications """
    print "making classification"

    # set up
    fileout = 'classification.csv'
    filelist = 'list_ran.csv'

    tab_header = at.Table.read(dir_batch+filelist, format='ascii.csv')
    tab_contam = classify_contamination(dir_batch, bandmag='r', min_modelMag=18)
    tab_detect = classify_detection(dir_batch)
    tab_resolv = classify_resolved(dir_batch)
    tab_bgtmask = classify_bgtmask(dir_batch)
    tab_patchy = classify_patchy(dir_batch)

    tab_class = at.hstack([tab_header, tab_contam, tab_detect, tab_resolv, tab_bgtmask, tab_patchy])

    tab_class.write(dir_batch+fileout, format='ascii.csv')


def classify_patchy(dir_batch):
    """
    classify things that are patchy. blob > blobctr
    """

    fin1 = dir_batch+'measureimg_iso_contours_blob.csv'
    fin2 = dir_batch+'measureimg_iso_contours_blobctr.csv'

    t1 = at.Table.read(fin1, format='ascii')
    t2 = at.Table.read(fin2, format='ascii')

    arr_patchy = t1['isoshape_area'] > t2['isoshape_area']

    tab_patchy = at.Table(data=[arr_patchy], names=['patchy', ], masked=False)
    return tab_patchy


def classify_bgtmask(dir_batch):
    """
    classify things that are greater than galaxy mask and psf mask 
    """
    # setup
    fin = dir_batch+'measureimg_iso_contours_blob_psfresid_galmasked_psfmasked.csv'
    tabin = at.Table.read(fin)

    arr_bgtmask = (tabin['isoshape_area'] > 0)
    tab_bgtmask = at.Table(data=[arr_bgtmask], names=['bgtmask', ], masked=False)
    return tab_bgtmask


def classify_resolved(dir_batch):
    """
    classify resolved 
    """
    # setup
    fin = dir_batch+'measureimg_iso_contours_blob_psfresid.csv'
    tabin = at.Table.read(fin)

    arr_resolved = (tabin['isoshape_area'] > 0)
    tab_resolv = at.Table(data=[arr_resolved], names=['resolved', ], masked=False)
    return tab_resolv


def classify_detection(dir_batch):
    """
    classify detection 

    Returns
    ----
    tab_detect: unmasked at table of one bool column:
        'detected'
    """
    # setup
    fin = dir_batch+'measureimg_iso_contours_blob.csv'
    tabin = at.Table.read(fin)

    arr_detected = (tabin['isoshape_area'] > 0)
    tab_detect = at.Table(data=[arr_detected], names=['detected', ], masked=False)
    return tab_detect


def classify_contamination(dir_batch, bandmag='r', min_modelMag=18):
    """
    classify contamination - whether there are close by object brighter than 
    'min_modelMag' at band 'bandmag'. 

    Params
    ------
    dir_batch, bandmag='r', min_modelMag=18

    Read Inputs
    ------
    dir_batch+'contaminants.csv'

    Returns
    ------
    tab_contam: masked at table w two columns: 
        'contamfree': True if no contamination 
        'contambygalonly': True if only contaminated by galaxy
        
    """
    # setup
    fin = dir_batch+'contaminants.csv'
    tabin = at.Table.read(fin)
    magtag = '_mag'+bandmag+str(min_modelMag)

    # calculate
    arr_contamfree = (tabin['n_contaminants'+magtag] == 0).data
    arr_contambygalonly = np.zeros(len(tabin), dtype=bool)
    for i, a in enumerate(tabin['types'+magtag]):
        if tabin['types'+magtag].mask[i]:  # masked
            arr_contambygalonly[i] = np.nan
        else: 
            arr_contambygalonly[i] = (set(str(a)) == {'3'}) # galaxy only

    # organize output
    tab_contam = at.Table(data=[arr_contamfree, arr_contambygalonly], names=['contamfree', 'contambygalonly'], masked=True)
    tab_contam['contamfree'].mask = tabin['n_contaminants'+magtag].mask
    tab_contam['contambygalonly'].mask = tabin['types'+magtag].mask

    return tab_contam


