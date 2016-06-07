# getzrange_batch.py
# ALS 2016/05/04

import numpy as np
from astropy.table import Table, vstack, hstack
import filtertools

def write_zranges(filenameout='', wline='OIII', nline='HaNIISII'):
    """
    write table zranges_batch.csv
    """
    list_kwargs=[{'bandline':'g','bandconti':'i','wline':wline,'nline':nline},
                {'bandline':'g','bandconti':'r','wline':wline,'nline':nline},
                {'bandline':'r','bandconti':'z','wline':wline,'nline':nline},
                {'bandline':'r','bandconti':'i','wline':wline,'nline':nline}]
    
    tabout=Table()

    for kwargs in list_kwargs:
        # calculation
        z0, z1=getzrange_batch(**kwargs)

        # set up table
        cols=['bandline','bandconti']
        tabheader=Table([[kwargs[col]] for col in cols],names=cols)
        tabdata=Table([[z0],[z1],[np.mean([z0,z1])]],names=['z0','z1','z_mean'])
        tabrow=hstack([tabheader,tabdata])
        if len(tabout)==0: tabout=tabrow
        else: tabout=vstack([tabout,tabrow])

    if filenameout=='':
        localpath=filtertools.getlocalpath()
        filenameout=localpath+'zranges_band_w'+wline+'_n'+nline+'.txt'
    tabout.write(filenameout,format='ascii.fixed_width',delimiter='')
    return tabout


def isinside(b, a0, a1):
    """
    return if b is inside the range between a0 and a1
    """
    a0,a1=np.sort([a0,a1])

    return ((b >= a0) and (b <= a1))


def getzrange_batch(bandline='r',bandconti='i',wline='OIII',nline='HaNIISII',wthreshold=0.6, nthreshold=0.2):
    """
    return the redshift range where
        bandline (r) contains all wline (OIII) and 
        contiband (i) contains non of nline (HaNIISII)

    """
    filename_wline='zrange_wline_'+wline+'_'+'%.1f'%wthreshold+'.txt'
    filename_nline='zrange_nline_'+nline+'_'+'%.1f'%nthreshold+'.txt'

    zranges_wline=filtertools.accessFile(filename=filename_wline)
    zranges_nline=filtertools.accessFile(filename=filename_nline)

    z0_bl, z1_bl = zranges_wline[zranges_wline['band']==bandline]['z0','z1'][0]
    z0_bc, z1_bc = zranges_nline[zranges_nline['band']==bandconti]['z0','z1'][0]

    bl_zmin, bl_zmax=np.sort([z0_bl, z1_bl])
    bc_zmin, bc_zmax=np.sort([z0_bc, z1_bc])


    bl_zmin, bl_zmax, bc_zmin, bc_zmax= max(0,bl_zmin), max(0,bl_zmax), max(0,bc_zmin), max(0,bc_zmax)

    bl_zmin_isinside=isinside(bl_zmin, bc_zmin, bc_zmax)
    bl_zmax_isinside=isinside(bl_zmax, bc_zmin, bc_zmax)

    if bl_zmin_isinside and not bl_zmax_isinside :
        z0=max(bl_zmin, bc_zmax)  
        z1=bl_zmax
    elif not bl_zmin_isinside and bl_zmax_isinside :
        z0=bl_zmin
        z1=min(bl_zmax, bc_zmin)
    else: 
        print "bl_zmin, bl_zmax, bc_zmin, bc_zmax:", str(bl_zmin), str(bl_zmax), str(bc_zmin), str(bc_zmax)
        raise ValueError("zrange errors -- scenario not considered")


    if z0>z1: raise ValueError("zrange errors")

    return z0, z1



# def getzrange_batch(bandline='r',bandconti='i',blue_edge='HaNIISII',red_edge='OIII'):
#     """
#     return the redshift range given the band set ups

#     blue_edge: which line determines the low z cut
#     red_edge: which line determines the high z cut
#     """
#     zrangesOIII=filtertools.accessFile(filename='OIIIredshiftrange0.6.txt')
#     zrangesHaNII=filtertools.accessFile(filename='HaNIIredshiftrange0.2.txt')
#     zrangesHaNIISII=filtertools.accessFile(filename='HaNIISIIredshiftrange0.2.txt')
#     # zrangesOIII=Table.read('OIIIredshiftrange0.6.txt',format='ascii')
#     # zrangesHaNII=Table.read('HaNIIredshiftrange0.2.txt',format='ascii')

#     # set low z cut
#     if blue_edge=='HaNII':
#         z1=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z2']
#     elif blue_edge=='HaNIISII':
#         z1=zrangesHaNIISII[zrangesHaNIISII['band']==bandconti][0]['z2']
#     elif blue_edge=='OIII':
#         z1=zrangesOIII[zrangesOIII['band']==bandline][0]['z1']
#     else: 
#         raise NameError('edge instruction not understood')

#     # set high z cut
#     if red_edge=='HaNII':
#         z2=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z1']
#     elif red_edge=='HaNIISII':
#         z2=zrangesHaNIISII[zrangesHaNIISII['band']==bandconti][0]['z1']
#     elif red_edge=='OIII':
#         z2=zrangesOIII[zrangesOIII['band']==bandline][0]['z2']
#     else: 
#         raise NameError('edge instruction not understood')

#     z1=max(0,z1)
#     z2=max(0,z2)

#     return z1, z2

# def write_zranges():
#     """
#     write table zranges_batch.csv
#     """

#     list_kwargs=[{'bandline':'g','bandconti':'i','blue_edge':'OIII','red_edge':'HaNIISII'},
#                 {'bandline':'g','bandconti':'r','blue_edge':'HaNIISII','red_edge':'OIII'},
#                 {'bandline':'r','bandconti':'z','blue_edge':'OIII','red_edge':'HaNIISII'},
#                 {'bandline':'r','bandconti':'i','blue_edge':'HaNIISII','red_edge':'OIII'}]
    
#     tabout=Table()

#     for kwargs in list_kwargs:
#         z0, z1=getzrange_batch(**kwargs)


#         # set up table
#         cols=['bandline','bandconti','blue_edge','red_edge']
#         tabheader=Table([[kwargs[col]] for col in cols],names=cols)

#         tabdata=Table([[z0],[z1],[np.mean([z0,z1])]],names=['z0','z1','z_mean'])

#         tabrow=hstack([tabheader,tabdata])
#         if len(tabout)==0: tabout=tabrow
#         else: tabout=vstack([tabout,tabrow])

#     localpath=filtertools.getlocalpath()
#     tabout.write(localpath+'zranges_band.txt',format='ascii.fixed_width',delimiter='')



#         z0, z1=getzrange_batch(**kwargs)
