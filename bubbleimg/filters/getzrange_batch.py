# getzrange_batch.py
# ALS 2016/05/04
import os
import numpy as np
import astropy.table as at

from . import filtertools


def write_zranges(filenameout='', wline='OIII', nline='HaNIISII', survey='sdss', wthreshold=0.6, nthreshold=0.2):
    """
    write table zranges_batch.csv
    """
    # set linetag
    linetag = '_w{wline}_n{nline}'.format(wline=wline, nline=nline)
    rangetag = '_w{wline}{wthreshold}_n{nline}{nthreshold}'.format(wline=wline, wthreshold=str(wthreshold), nline=nline, nthreshold=str(nthreshold))

    # set output file name
    if filenameout == '':
        localpath = filtertools.getlocalpath()
        filepathout = os.path.join(localpath+survey, 'zranges_band{}.txt'.format(rangetag))
    else: 
        filepathout = filenameout 

    # read in bandconfig
    tab_bandconfigs = filtertools.accessFile(filename='bandconfigs'+linetag+'.txt', survey=survey)
    tab_bandconfigs.meta['comments']=[]

    # make table
    tabout = at.Table()
    for tbandconfig in tab_bandconfigs:
        dictbandconfig = {key: tbandconfig[key][0] for key in tbandconfig.columns}
        dictbandconfig.update({'wline':wline,'nline':nline})
        # calculation
        z0, z1 = getzrange_batch(survey=survey, wthreshold=wthreshold, nthreshold=nthreshold, **dictbandconfig)
        # set up table

        if type(z0) is not list:
            tabdata = at.Table([[z0], [z1], [np.mean([z0, z1])]],names=['z0', 'z1', 'z_mean'])
            tabout = append_tabdata_to_tabout(tabdata, tbandconfig, tabout)            
        else:
            for i in range(len(z0)):
                tabdata = at.Table([[z0[i]], [z1[i]], [np.mean([z0[i], z1[i]])]],names=['z0', 'z1', 'z_mean'])
                tabout = append_tabdata_to_tabout(tabdata, tbandconfig, tabout)

    # output
    tabout.write(filepathout,format='ascii.fixed_width',delimiter='', overwrite=True)
    return tabout


def append_tabdata_to_tabout(tabdata, tbandconfig, tabout):
    tabrow = at.hstack([tbandconfig, tabdata])

    if len(tabout)==0: 
        tabout=tabrow
    else: 
        tabout=at.vstack([tabout,tabrow])
    return tabout


def getzrange_batch(bandline='r', bandconti='i', wline='OIII', nline='HaNIISII', wthreshold=0.6, nthreshold=0.2, survey='sdss'):
    """
    return the redshift range where
        bandline (r) contains all wline (OIII) and 
        contiband (i) contains non of nline (HaNIISII)

    Return
    ------
    z0 (float or list of floats)
    z1 (float or list of floats)
    """
    print(("calculating zrange_batch for "+bandline+'-'+bandconti))
    filename_wline = 'zrange_wline_'+wline+'_'+'%.1f'%wthreshold+'.txt'
    filename_nline = 'zrange_nline_'+nline+'_'+'%.1f'%nthreshold+'.txt'

    zranges_wline = filtertools.accessFile(filename=filename_wline, survey=survey)
    zranges_nline = filtertools.accessFile(filename=filename_nline, survey=survey)

    z0_bl, z1_bl = zranges_wline[zranges_wline['band']==bandline]['z0','z1'][0]
    z0_bc, z1_bc = zranges_nline[zranges_nline['band']==bandconti]['z0','z1'][0]

    bl_zmin, bl_zmax = np.sort([z0_bl, z1_bl])
    bc_zmin, bc_zmax = np.sort([z0_bc, z1_bc])


    bl_zmin, bl_zmax, bc_zmin, bc_zmax = max(0,bl_zmin), max(0,bl_zmax), max(0,bc_zmin), max(0,bc_zmax)

    z0, z1 = determine_non_excluded_zrange(in_zmin=bl_zmin, in_zmax=bl_zmax, ex_zmin=bc_zmin, ex_zmax=bc_zmax)

    return z0, z1


def determine_non_excluded_zrange(in_zmin, in_zmax, ex_zmin, ex_zmax):
    """

    return z range (z0, z1) that is within (in_zmin, in_zmax) and outside (ex_zmin, ex_zmax). 
    in case where there is more than 1 range that satisfy the criteria then return list ((z0_1, z0_2), (z1_1, z1_2))
    """
    in_zmin, in_zmax = np.sort([in_zmin, in_zmax])
    ex_zmin, ex_zmax = np.sort([ex_zmin, ex_zmax])

    in_zmin_isinside = isinside(in_zmin, ex_zmin, ex_zmax)
    in_zmax_isinside = isinside(in_zmax, ex_zmin, ex_zmax)

    # there is no overlap between the two ranges (in_zmin, in_zmax) and (ex)
    isnooverlap = ((ex_zmin < in_zmin) & (ex_zmax < in_zmin)) | ((ex_zmin > in_zmax) & (ex_zmax > in_zmax))

    if in_zmin_isinside and not in_zmax_isinside: # lower end of inz is excluded
        z0 = max(in_zmin, ex_zmax)  
        z1 = in_zmax

    elif not in_zmin_isinside and in_zmax_isinside: # higher end of inz is excluded
        z0 = in_zmin
        z1 = min(in_zmax, ex_zmin)

    elif (not in_zmin_isinside) and (not in_zmax_isinside) and (not isnooverlap): # middle part of inz is excluded
        z0 = [in_zmin, ex_zmax]
        z1 = [ex_zmin, in_zmax]

    elif in_zmin_isinside and in_zmax_isinside and (not isnooverlap): # entirely excluded
        z0 = np.nan
        z1 = np.nan

    elif isnooverlap:  # none of inz is excluded
        z0 = in_zmin
        z1 = in_zmax

    else: # should not happen
        print(("in_zmin, in_zmax, ex_zmin, ex_zmax:", str(in_zmin), str(in_zmax), str(ex_zmin), str(ex_zmax)))
        raise ValueError("zrange errors -- scenario not considered")

    if type(z0) is not list:
        if z0 > z1: 
            raise ValueError("[filters.getzrange_batch] zrange errors")
    elif type(z0) is list:
        if len(z0) != len(z1):
            raise ValueError("[filters.getzrange_batch] different z size")
        for i in range(len(z0)):
            if z0[i] > z1[i]:
                raise ValueError("[filters.getzrange_batch] zrange errors")

    return z0, z1


def isinside(b, a0, a1):
    """
    return if b is inside the range between a0 and a1
    """
    a0,a1=np.sort([a0,a1])

    return ((b >= a0) and (b <= a1))




