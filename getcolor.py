# getcolor.py
# ALS 2015/05/07
"""
PURPOSE: record and plot SDSS colors on the SDSS stamp images
"""

from pylab import *
import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/')
import galaxy


from scipy import interpolate
from astropy.table import Table
from astropy.io import fits
import astropy.units as u
# import class_obsobj
# reload(class_obsobj)
# from class_obsobj import obsobj

dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/Magellan/'
filelist='list_Magellan.txt'

nanomaggy = u.def_unit('nanomaggy', 3.631e-6*u.Jy)
u.add_enabled_units([nanomaggy])
u.nanomaggy=nanomaggy


def poswSlitFluxDensity(pos,bands=['u','g','r','i','z'],imagepixsize=0.396,slitpixsize=0.6,slitlength=20.,toplot=True):
	"""
	PURPOSE: obtain color along magellan slit position given "pos" object

	INPUT: 
		pos: position object
		imagepixsize: pix size of image [arcsec]
		slitpixsize:  pix size of slit  [arcsec]
		slitlength:   length of slit    [arcsec]

	OUTPUT: 
		bintable: 

		binindex  ycoord    Fnu_g_nmaggy      Fnu_r_nmaggy    Fnu_i_nmaggy  
		          arcsec                                                    
		 int64   float64      float64           float64         float64     
		-------- ------- ------------------ --------------- ----------------
	"""

	# setup
	dir_obj=dir_data+'M'+pos.OBJID+'/'
	fileout=dir_obj+'Obj'+pos.OBJID+'P'+pos.POS+'_slitfnu'+'_bsw'+str(slitpixsize)+'ars'

	# initialize bintable, containing columns ("binindex" [pix]  "ycoord" [arcsec])
	nbin=int(floor(slitlength/slitpixsize / 2.) * 2 +1) # number of bins rounded to closest odd number
	binindexes=arange(nbin)-(nbin-1)/2
	binycoord=binindexes*slitpixsize
	bintable=Table(data=[binindexes,binycoord],names=['binindex','ycoord'])
	bintable['ycoord'].unit=u.arcsec

	# writing flux density
	for band in bands:	# band='r'
		# create row in bintable to store flux density
		bintable['Fnu_'+band+'_nmaggy']=0.#*u.nanomaggy

		# readin image
		filename=dir_obj+'stamp-'+band+'.fits'
		hdulist=fits.open(filename)
		image=hdulist[0].data
		header=hdulist[0].header
		# image spec
		nx,ny=image.shape
		x0,y0=header['CRPIX1'],header['CRPIX2']
		# spline 
		f=interpolate.interp2d(arange(nx), arange(ny), image, kind='cubic')
		# relative PA (alpha)
		imagePA=getimagePA_from_header(hdulist[0].header)
		alpha=pos.PA-imagePA

		for i in range(nbin):
			r=bintable['ycoord'][i]
			x=-r*sin(radians(alpha))/imagepixsize+x0 # (x,y) vector from image center to the bin [image pix]
			y= r*cos(radians(alpha))/imagepixsize+y0

			bintable['Fnu_'+band+'_nmaggy'][i]=f(x,y)
			# print x, y

	# output
	bintable.write(fileout+'.fits',format='fits',overwrite=True)
	bintable.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')

	if toplot:
		clf()
		for band in bands:	
			plot(bintable['ycoord'],bintable['Fnu_'+band+'_nmaggy'],label=band)
		axhline(0,color='black')
		legend()
		xlabel('ycoord [arcsec]')
		ylabel('Fnu [nanomaggy]')
		savefig(fileout+'.pdf')



def poswSlitColor(pos,bands=['u','g','r','i','z'],imagepixsize=0.396,slitpixsize=0.6,slitlength=20.,toplot=True):
	"""
	PURPOSE: digest slitfluxdensity table and write slitcolor table that contains color, intensity, etc. 

	OUTPUS: 
		bintable: 

		binindex  ycoord    Fnu_g_nmaggy      Fnu_r_nmaggy    Fnu_i_nmaggy         Fnu_g              Fnu_r       ...       Snu_r             Snu_i          g-r_mag         r-i_mag         r_over_g      r_over_i   
		          arcsec                                                      erg / (cm2 Hz s)   erg / (cm2 Hz s) ... erg / (Hz kpc2 s) erg / (Hz kpc2 s)      mag             mag                                    
		 int64   float64      float64           float64         float64           float64            float64      ...      float64           float64         float64         float64         float64       float64    
		-------- ------- ------------------ --------------- ---------------- ------------------ ----------------- ... ----------------- ----------------- -------------- ---------------- ------------- --------------

	"""
	# setup
	dir_obj=dir_data+'M'+pos.OBJID+'/'
	filein=dir_obj+'Obj'+pos.OBJID+'P'+pos.POS+'_slitfnu'+'_bsw'+str(slitpixsize)+'ars.fits'
	fileout=dir_obj+'Obj'+pos.OBJID+'P'+pos.POS+'_slitcolor'+'_bsw'+str(slitpixsize)+'ars'

	z=pos.outer.outer.ZSTAR
	bandratiopairs=[['r','u'],['r','g'],['r','i'],['r','z']]
	extcolorpairs=[['u','r'],['r','z']]
	noisethreshold=0.02 # nanomaggy
	seterr(all='ignore')

	# read in
	bintable=Table.read(filein,format='fits')

	#===== operation
	# define noisy bins

	selnoisy=zeros(len(bintable))
	for band in bands:	
		selnoisy=logical_or(selnoisy, bintable['Fnu_'+band+'_nmaggy']<noisethreshold)

	# create columns Fnu in units of [erg / (cm2 Hz s)]
	for band in bands:	
		bintable['Fnu_'+band]=(bintable['Fnu_'+band+'_nmaggy']*u.nanomaggy).to(u.erg/u.s/u.cm**2/u.Hz)

	# create columns Snu in units of [erg / (cm2 Hz s)]
	for band in bands:	
		# WARNING this conversion has to be checked
		bintable['Snu_'+band]=(bintable['Fnu_'+band]*4.*pi*u.sr/(imagepixsize*u.arcsec)**2*(1.+z)**3).to(u.erg/u.s/u.kpc**2/u.Hz)

	# creat columns for colors
	for ib in range(len(bands)-1):
		bintable[bands[ib]+'-'+bands[ib+1]+'_mag']=u.mag*(-2.5*log10(bintable['Fnu_'+bands[ib]+'_nmaggy']/bintable['Fnu_'+bands[ib+1]+'_nmaggy']))
		bintable[bands[ib]+'-'+bands[ib+1]+'_mag'][selnoisy]=nan

	for bandpair in extcolorpairs:
		bintable[bandpair[0]+'-'+bandpair[1]+'_mag']=u.mag*(-2.5*log10(bintable['Fnu_'+bandpair[0]+'_nmaggy']/bintable['Fnu_'+bandpair[1]+'_nmaggy']))
		bintable[bandpair[0]+'-'+bandpair[1]+'_mag'][selnoisy]=nan


	# creat columns for ratio
	for bandpair in bandratiopairs:
		bintable[bandpair[0]+'_over_'+bandpair[1]]=bintable['Fnu_'+bandpair[0]+'_nmaggy']/bintable['Fnu_'+bandpair[1]+'_nmaggy']
		bintable[bandpair[0]+'_over_'+bandpair[1]][selnoisy]=nan



	# output
	bintable.write(fileout+'.fits',format='fits',overwrite=True)
	bintable.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')

	if toplot:
		clf()
		for ib in range(len(bands)-1):
			plot(bintable['ycoord'],bintable[bands[ib]+'-'+bands[ib+1]+'_mag'],marker='.',label=bands[ib]+'-'+bands[ib+1])
		legend()
		xlabel('ycoord [arcsec]')
		ylabel('color [mag]')
		savefig(fileout+'_colors.pdf')

		clf()
		for bandpair in bandratiopairs:
			plot(bintable['ycoord'],bintable[bandpair[0]+'_over_'+bandpair[1]],marker='.',label=bandpair[0]+'/'+bandpair[1])
		legend()
		xlabel('ycoord [arcsec]')
		ylabel('ratio')
		ylim(0.,4.)
		savefig(fileout+'.pdf')


def plotallcolors():
	"""
	PURPOSE: plot all ratios of all slits
	"""
	# setup
	troadmap=galaxy.gettroadmap()
	troadmap.sort('OBJID')
	bandratiopairs=[['r','g'],['r','i']]


	for bandpair in bandratiopairs:

		clf()
		for i in range(len(troadmap)):
			OBJID,POS=troadmap['OBJID','POS'][i]
			OBJID=OBJID.astype('str')
			POS=POS.astype('str')
			filein=dir_data+'M'+OBJID+'/'+'Obj'+OBJID+'P'+POS+'_slitcolor'+'_bsw'+str(slitpixsize)+'ars.fits'

			# readin
			bintable=Table.read(filein,format='fits')
			# plot
			plot(bintable['ycoord'],bintable[bandpair[0]+'_over_'+bandpair[1]],marker='.',label=OBJID)

		legend(prop={'size':8})
		xlabel('ycoord [arcsec]')
		ylabel('ratio '+bandpair[0]+'/'+bandpair[1])
		savefig(dir_data+'ratio_'+bandpair[0]+'_over_'+bandpair[1]+'.pdf')


def getimagePA_from_header(header):
	"""
	PURPOSE: getting the PA (y-axis east of north) of the image from image header. 

	NOTES:  sintheta = dDec/dx = CD21/sqrt(CD21^2+CD11^2)
		   -costheta = dRA/dx  = CD11/sqrt(CD21^2+CD11^2)

		   	costheta = dDec/dy = CD22/sqrt(CD22^2+CD12^2)
		   	sintheta = dRA/dy  = CD12/sqrt(CD22^2+CD12^2)
	"""
	sintheta=header['CD2_1']/sqrt(header['CD2_1']**2+header['CD1_1']**2)

	if header['CD1_1']<=0:
		pa=degrees(arcsin(sintheta))
	else: 
		pa=180.-degrees(arcsin(sintheta))

	return pa




def test():
	
	OBJID='2100'
	images=np.load('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/Magellan/M'+OBJID+'/images_stamp.npy')
	# d=ds9.ds9('8070198b:49216')


	d.set_np2arr(images[1]/sqrt(images[0]*images[2]))

	a=images[1]/images[0]#sqrt(images[0]*images[2])

	a[a>5.]=0.
	a[a<1]=0.
	a[~isfinite(a)]=0.

	d.set_np2arr(a)
