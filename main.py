# main.py
# ALS 2015/05/07
"""
automatically make the stamp images, color images for Magellan targets
"""

import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/')
import galaxy


import makeimages
reload(makeimages)
import getcolor
reload(getcolor)

def main():

	makeimages.makeStampsImages()

	makeimages.makeHumVIgriimages()



	galaxy.posIterator(getcolor.poswSlitFluxDensity)
	galaxy.posIterator(getcolor.poswSlitColor)
	
if __name__ == '__main__':
	main()

