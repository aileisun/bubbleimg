# Magellan_main.py
# ALS reorganized 2015/08/10

"""
 For Magellan targets automatically make the stamp images, color images, and extract color along slits

"""


from pylab import *
from astropy.table import Table, hstack, vstack
import os

import class_obsobj
reload(class_obsobj)
from class_obsobj import obsobj

import alignstamp
reload(alignstamp)

import subtractimg
reload(subtractimg)

import fromspec
reload(fromspec)


import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/')
import galaxy

import magellan_getslit
reload(magellan_getslit)

import imagedisp_util
reload(imagedisp_util)

import measurenebula
reload(measurenebula)


dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/magellan/'



def main():
	"""
	automatically make the stamp images, color images, and extract color along slits for Magellan targets

	products saved to self.dir_obj=magellan_dir_data+'M'+str(self.galaxy.OBJID)+'/'
	"""

	#==== set up
	filelist=dir_data+'list_Magellan.txt'
	listmagellan=Table.read(filelist,format='ascii')

	for i in range(len(listmagellan)):
		# declare obj
		obj=obsobj(listmagellan[i],catalog='magellan')
		print obj.galaxy.OBJID

		# make stamps		
		alignstamp.objw_makeall(obj)
		# make color images
		imagedisp_util.objw_HumVIgriimages(obj)
		# subtract stamps and make OIII maps
		subtractimg.objw_makeall(obj)

	# Extract info from stamps along Magellan slits
	galaxy.posIterator(magellan_getslit.poswSlitFluxDensity)

	kwargs={'linetag':'lOIII5008'}
	galaxy.posIterator(magellan_getslit.poswSlitLineIntensity,**kwargs)

	galaxy.posIterator(magellan_getslit.poswSlitAll)

	measurenebula.write_measureISO_lOIII5008(dir_data=dir_data,catalog='magellan',filelist='list_Magellan.txt')

	

if __name__ == '__main__':
	main()
