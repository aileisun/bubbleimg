# getzrange_batch.py
# ALS 2016/05/04

import numpy as np
import astropy.table as at
import filtertools

def write_zranges(filenameout='', wline='OIII', nline='HaNIISII', survey='sdss'):
    """
    write table zranges_batch.csv
    """
    # set linetag
    linetag = '_w'+wline+'_n'+nline

    # set output file name
    if filenameout=='':
        localpath = filtertools.getlocalpath()
        filepathout = localpath+survey+'/'+'zranges_band'+linetag+'.txt'
    else: 
        filepathout = filenameout 

    # read in bandconfig
    tab_bandconfigs = filtertools.accessFile(filename='bandconfigs'+linetag+'.txt', survey=survey)
    tab_bandconfigs.meta['comments']=[]

    # make table
    tabout=at.Table()
    for tbandconfig in tab_bandconfigs:
        dictbandconfig = {key: tbandconfig[key][0] for key in tbandconfig.columns}
        dictbandconfig.update({'wline':wline,'nline':nline})
        # calculation
        z0, z1=getzrange_batch(survey=survey, **dictbandconfig)
        # set up table
        cols=['bandline','bandconti']
        tabdata = at.Table([[z0],[z1],[np.mean([z0,z1])]],names=['z0','z1','z_mean'])
        tabrow = at.hstack([tbandconfig,tabdata])
        if len(tabout)==0: tabout=tabrow
        else: tabout=at.vstack([tabout,tabrow])

    # output
    tabout.write(filepathout,format='ascii.fixed_width',delimiter='', overwrite=True)
    return tabout


def getzrange_batch(bandline='r',bandconti='i',wline='OIII',nline='HaNIISII',wthreshold=0.6, nthreshold=0.2, survey='sdss'):
    """
    return the redshift range where
        bandline (r) contains all wline (OIII) and 
        contiband (i) contains non of nline (HaNIISII)

    """
    filename_wline = 'zrange_wline_'+wline+'_'+'%.1f'%wthreshold+'.txt'
    filename_nline = 'zrange_nline_'+nline+'_'+'%.1f'%nthreshold+'.txt'

    zranges_wline = filtertools.accessFile(filename=filename_wline, survey=survey)
    zranges_nline = filtertools.accessFile(filename=filename_nline, survey=survey)

    z0_bl, z1_bl = zranges_wline[zranges_wline['band']==bandline]['z0','z1'][0]
    z0_bc, z1_bc = zranges_nline[zranges_nline['band']==bandconti]['z0','z1'][0]

    bl_zmin, bl_zmax=np.sort([z0_bl, z1_bl])
    bc_zmin, bc_zmax=np.sort([z0_bc, z1_bc])


    bl_zmin, bl_zmax, bc_zmin, bc_zmax= max(0,bl_zmin), max(0,bl_zmax), max(0,bc_zmin), max(0,bc_zmax)

    bl_zmin_isinside=isinside(bl_zmin, bc_zmin, bc_zmax)
    bl_zmax_isinside=isinside(bl_zmax, bc_zmin, bc_zmax)

    isnooverlap = ((bc_zmin < bl_zmin) & (bc_zmax < bl_zmin)) | ((bc_zmin > bl_zmax) & (bc_zmax > bl_zmax))

    if bl_zmin_isinside and not bl_zmax_isinside :
        z0=max(bl_zmin, bc_zmax)  
        z1=bl_zmax
    elif not bl_zmin_isinside and bl_zmax_isinside :
        z0=bl_zmin
        z1=min(bl_zmax, bc_zmin)
    elif isnooverlap: 
        z0=bl_zmin
        z1=bl_zmax
    else: 
        print "bl_zmin, bl_zmax, bc_zmin, bc_zmax:", str(bl_zmin), str(bl_zmax), str(bc_zmin), str(bc_zmax)
        raise ValueError("zrange errors -- scenario not considered")


    if z0>z1: raise ValueError("zrange errors")

    return z0, z1


def isinside(b, a0, a1):
    """
    return if b is inside the range between a0 and a1
    """
    a0,a1=np.sort([a0,a1])

    return ((b >= a0) and (b <= a1))






    # bandconfigs = 

    # # define operations
    # if survey == 'sdss': 
    #     bandconfigs=[{'bandline':'g','bandconti':'i'},
    #                 {'bandline':'g','bandconti':'r'},
    #                 {'bandline':'r','bandconti':'z'},
    #                 {'bandline':'r','bandconti':'i'},
    #                 {'bandline':'i','bandconti':'z'},
    #                 ]
    # elif survey == 'hsc': 
    #     bandconfigs=[{'bandline':'g','bandconti':'i'},
    #                 {'bandline':'g','bandconti':'r'},
    #                 {'bandline':'r','bandconti':'z'},
    #                 {'bandline':'r','bandconti':'i'},
    #                 {'bandline':'i','bandconti':'z'},
    #                 {'bandline':'z','bandconti':'y'},
    #                 ]
    # elif survey == 'cfht-hsc-ukirt': 
    #     bandconfigs=[
    #                 {'bandline':'g','bandconti':'i'},
    #                 {'bandline':'g','bandconti':'r'},
    #                 {'bandline':'r','bandconti':'z'},
    #                 {'bandline':'r','bandconti':'y'},
    #                 {'bandline':'r','bandconti':'i'},
    #                 {'bandline':'i','bandconti':'y'},
    #                 {'bandline':'i','bandconti':'z'},
    #                 {'bandline':'z','bandconti':'j'},
    #                 {'bandline':'j','bandconti':'k'},
    #                 {'bandline':'k','bandconti':'j'},
    #                 ]
    # else: 
    #     raise NameError('survey name not recognized. ')
    