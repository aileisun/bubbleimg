# filterzrange.py
# ALS 2016/05/02
"""
Determine the redshift range of a given line in a given filter. The outputs 
are .txt files: 

	filterboundary_0.2.txt
	filterboundary_0.6.txt
	filterboundary_0.8.txt
	HaNIIredshiftrange0.2.txt
	OIIIredshiftrange0.6.txt

"""
# from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Row
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy.interpolate import UnivariateSpline

import filtertools
import getzrange_band

def main():
	"""
	make files: 
		filterboundary_0.2.txt
		filterboundary_0.6.txt
		filterboundary_0.8.txt
		HaNIIredshiftrange0.2.txt
		OIIIredshiftrange0.6.txt
	"""
	for threshold in [0.2, 0.6, 0.8]:
		findFilterBounday(threshold=threshold,toplot=False)

	# write files
	findRedshiftRange_OIIIs(threshold=0.6)
	findRedshiftRange_HaNII(threshold=0.2)
	findRedshiftRange_HaNIISII(threshold=0.2)
	
	getzrange_band.write_zranges()
		

def findFilterBounday(threshold=0.6,toplot=True):
	"""
	PURPOSE: write file filterboundary.txt to record filter boundary defined by 80% of maximum throughput
	PARAMETERS:
			threshold=0.6 (float)
			toplot=True
	"""

	# fileout
	fileout='filterboundary_'+'%.1f'%threshold+'.txt'

	tabout=Table([[],[],[],],names=('band','w1','w2'),dtype=('string','int','int'))

	for band in ['u','g','r','i','z']:
		# readin filter function
		spec,lcoord=filtertools.getFilterResponseFunc(band=band)


		spl=UnivariateSpline(lcoord, spec/max(spec)-threshold,s=0)
		roots=spl.roots()
		# print roots

		# solve for root
		f = interp1d(lcoord, spec/max(spec),kind='linear',bounds_error=False,fill_value=0.)
		tosolve=lambda x: f(x)-threshold

		x1=fsolve(tosolve,x0=roots[0])
		x2=fsolve(tosolve,x0=roots[-1])

		tabout.add_row([band,x1,x2])

		if toplot:
			plt.clf()
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.plot(lcoord,spec/max(spec),color='black')
			ax.axhline(threshold)
			ax.axvline(x1)
			ax.axvline(x2)
			ax.set_xlim(min(lcoord),max(lcoord))
			# raw_input("Enter...")
			fig.savefig(band+'_'+'%.1f'%threshold+'.pdf')

	tabout.write(fileout,format='ascii.fixed_width',delimiter='')



def findRedshiftRange_OIIIs(threshold=0.6):
	"""
	Find the right red shift range such that both [OIII] 4960 and 5008 are within 
	throughput > 0.6*max range of r band. 
	"""
	fileprefix='OIII'
	l1=getllambda(ion='OIII',lid=4960,vacuum=True) 
	l2=getllambda(ion='OIII',lid=5008,vacuum=True) 

	lmin=min(l1,l2)
	lmax=max(l1,l2)

	findzrange_line(fileprefix,lmin,lmax,threshold=threshold)



def findRedshiftRange_HaNII(threshold=0.6):
	"""
	Find the right red shift range such that both [OIII] 4960 and 5008 are within 
	throughput > 0.6*max range of r band. 
	"""
	fileprefix='HaNII'
	l1=getllambda(ion='Ha')[0]
	l2,l3=getllambda(ion='NII')

	lmin=min(l1,l2,l3)
	lmax=max(l1,l2,l3)

	findzrange_line(fileprefix,lmin,lmax,threshold=threshold)



def findRedshiftRange_HaNIISII(threshold=0.6):
	"""
	Find the right red shift range such that both [OIII] 4960 and 5008 are within 
	throughput > 0.6*max range of r band. 
	"""
	fileprefix='HaNIISII'
	l1=getllambda(ion='Ha')[0]
	l2,l3=getllambda(ion='NII')
	l4,l5=getllambda(ion='SII')[1:]

	lmin=min(l1, l2, l3, l4, l5)
	lmax=max(l1, l2, l3, l4, l5)

	findzrange_line(fileprefix,lmin,lmax,threshold=threshold)




def findzrange_line(fileprefix,lmin,lmax,threshold=0.6):
	"""
	Find the right red shift range such that Halpha and two NII lines are within
	throughput > 0.6*max range of r band. 
	"""

	# setups
	fileout=fileprefix+'redshiftrange'+'%.1f'%threshold+'.txt'
	tabout=Table([[],[],[],],names=('band','z1','z2'),dtype=('string','float','float'))

	
	for band in ['u','g','r','i','z']:

		# get filter boundary
		fb=Table.read('filterboundary_'+str(threshold)+'.txt',format='ascii')
		w1,w2=fb[fb['band']==band]['w1','w2'][0]

		# calculate corresponding redshift
		z1=w1/lmin-1
		z2=w2/lmax-1
		tabout.add_row([band,z1,z2])

	# output
	tabout.write(fileout,format='ascii.fixed_width',delimiter='')
	return tabout



def getllambda(ion='OIII',lid=0,vacuum=True):
	"""
	PURPOSE: Get the vaccum or air wavelength of a line in Angstrom. 

	INPUT:  getllambda(ion='OIII',lid=0,vacuum=True)
		ion: e.g., ion='OIII'
		lid:  e.g., lid=5008. If lid not specified then all lines of that ion will be returned. 
		vacuum: True if want vacuum wavelength, false if want air. 

	OUTPUT: line wavelengths in Angstrom
	"""
	lid=int(lid)
	# e.g. ion='OIII', lid=5008
	linelist=Table.read('linelist.txt',format='ascii',delimiter='\t')
	if lid ==0:
		lam=linelist[linelist['Line']==ion]['Wavelength'].data
	else:
		if (linelist['Line']==ion).sum()==1:
			lam=linelist[linelist['Line']==ion]['Wavelength'].data
		elif (linelist['Line']==ion).sum()>1:
			lam=linelist[np.all([linelist['Line']==ion,linelist['Identifier']==lid],axis=0)]['Wavelength'].data
		else:
			raise NameError("[llambda] no line found")

	if vacuum:
		return lam
	else:
		from PyAstronomy import pyasl
		return pyasl.vactoair2(lam)


if __name__ == '__main__':
	main()

