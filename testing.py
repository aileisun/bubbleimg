test.py


from astropy.table import Table

import run_batch

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/')
import catalogue_util


list_torun=Table()

list_torun['RA']=[35.11413]
list_torun['DEC']=[1.23540]
list_torun['SDSSNAME']=['SDSS J022027.39+011407.4']
list_torun['Z']=[0.290]



dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/test/'

bandline='r'
bandconti='i'


run_batch.run_batch(list_torun,dir_data, bandline='r',bandconti='z',batch='',catalog='mullaney',tojoinmullaney=True,suffix='')


filemullaney='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/catalogue/ALPAKA_extended.fits'


tmullaney=Table.read(filemullaney,format='fits')

catalogue_util.selectRADEC(tmullaney,list_torun)
# NO MATCH 35.11413 1.2354
# This obj is not in Mullaney


