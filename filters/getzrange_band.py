# getzrange_band.py
# ALS 2016/05/04

import numpy as np
from astropy.table import Table, vstack, hstack
import filtertools

def write_zranges():
    """
    write table zranges_batch.csv
    """

    list_kwargs=[{'bandline':'g','bandconti':'i','blue_edge':'OIII','red_edge':'HaNIISII'},
                {'bandline':'g','bandconti':'r','blue_edge':'HaNIISII','red_edge':'OIII'},
                {'bandline':'r','bandconti':'z','blue_edge':'OIII','red_edge':'HaNIISII'},
                {'bandline':'r','bandconti':'i','blue_edge':'HaNIISII','red_edge':'OIII'}]
    
    tabout=Table()

    for kwargs in list_kwargs:
        z0, z1=getzrange(**kwargs)


        # set up table
        cols=['bandline','bandconti','blue_edge','red_edge']
        tabheader=Table([[kwargs[col]] for col in cols],names=cols)

        tabdata=Table([[z0],[z1],[np.mean([z0,z1])]],names=['z0','z1','z_mean'])

        tabrow=hstack([tabheader,tabdata])
        if len(tabout)==0: tabout=tabrow
        else: tabout=vstack([tabout,tabrow])

    localpath=filtertools.getlocalpath()
    tabout.write(localpath+'zranges_band.txt',format='ascii.fixed_width',delimiter='')



def getzrange(bandline='r',bandconti='i',blue_edge='HaNIISII',red_edge='OIII'):
    """
    return the redshift range given the band set ups

    blue_edge: which line determines the low z cut
    red_edge: which line determines the high z cut
    """
    zrangesOIII=filtertools.accessFile(filename='OIIIredshiftrange0.6.txt')
    zrangesHaNII=filtertools.accessFile(filename='HaNIIredshiftrange0.2.txt')
    zrangesHaNIISII=filtertools.accessFile(filename='HaNIISIIredshiftrange0.2.txt')
    # zrangesOIII=Table.read('OIIIredshiftrange0.6.txt',format='ascii')
    # zrangesHaNII=Table.read('HaNIIredshiftrange0.2.txt',format='ascii')

    # set low z cut
    if blue_edge=='HaNII':
        z1=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z2']
    elif blue_edge=='HaNIISII':
        z1=zrangesHaNIISII[zrangesHaNIISII['band']==bandconti][0]['z2']
    elif blue_edge=='OIII':
        z1=zrangesOIII[zrangesOIII['band']==bandline][0]['z1']
    else: 
        raise NameError('edge instruction not understood')

    # set high z cut
    if red_edge=='HaNII':
        z2=zrangesHaNII[zrangesHaNII['band']==bandconti][0]['z1']
    elif red_edge=='HaNIISII':
        z2=zrangesHaNIISII[zrangesHaNIISII['band']==bandconti][0]['z1']
    elif red_edge=='OIII':
        z2=zrangesOIII[zrangesOIII['band']==bandline][0]['z2']
    else: 
        raise NameError('edge instruction not understood')

    z1=max(0,z1)
    z2=max(0,z2)

    return z1, z2
