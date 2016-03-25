# getslit_Magellan.py
# ALS 2015/05/07
"""
PURPOSE: put artificial slits on the SDSS images that represent Magellan observations, extract broad-band-constructed line intensities along those 
	slits, and compare with the Magellan observations. 

"""

from pylab import *
import os
import sys
sys.path.append('/Users/aisun/Documents/Astro/Thesis/followups/Magellan/2014June/data_v2/analysis/')
import galaxy


from scipy import interpolate
from astropy.table import Table, join
from astropy.io import fits
import astropy.units as u
# import class_obsobj
# reload(class_obsobj)
# from class_obsobj import obsobj

dir_data='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/magellan/'

nanomaggy = u.def_unit('nanomaggy', 3.631e-6*u.Jy)
u.add_enabled_units([nanomaggy])
u.nanomaggy=nanomaggy



# def test():
	
# 	OBJID='2100'
# 	images=np.load('/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/data/Magellan/M'+OBJID+'/images_stamp.npy')
# 	# d=ds9.ds9('8070198b:49216')


# 	d.set_np2arr(images[1]/sqrt(images[0]*images[2]))

# 	a=images[1]/images[0]#sqrt(images[0]*images[2])

# 	a[a>5.]=0.
# 	a[a<1]=0.
# 	a[~isfinite(a)]=0.

# 	d.set_np2arr(a)



def poswSlitAll(pos,bands=['u','g','r','i','z'],imagepixsize=0.396,slitpixsize=0.6,slitlength=20.,toplot=True):
	"""
	PURPOSE: 
		summarize slitfluxdensity and slit line intensity table and write slitall table 
		that contains all info like color, intensity, etc. 

	WRITE OUTPUT:  
		e.g., Obj2100P1_slitcolor_bsw0.6ars.fits

		content: bintable

		binindex  ycoord      I_lOIII5008         Fnu_g_nmaggy      Fnu_r_nmaggy    Fnu_i_nmaggy         Fnu_g              Fnu_r       ...       Snu_r             Snu_i          g-r_mag         r-i_mag         r_over_g      r_over_i   
		          arcsec  erg / (arcsec2 cm2 s)                                                      erg / (cm2 Hz s)   erg / (cm2 Hz s) ... erg / (Hz kpc2 s) erg / (Hz kpc2 s)      mag             mag                                    
		 int64   float64         float64             float64           float64         float64           float64            float64      ...      float64           float64         float64         float64         float64       float64    
		-------- -------  --------------------- ------------------ --------------- ---------------- ------------------ ----------------- ... ----------------- ----------------- -------------- ---------------- ------------- --------------

	"""
	# setup
	dir_obj=dir_data+'M'+pos.OBJID+'/'
	# output
	fileout=dir_obj+'Obj'+pos.OBJID+'P'+pos.POS+'_slitall'+'_bsw'+str(slitpixsize)+'ars'
	# input
	linetag='lOIII5008'
	filein_lI=dir_obj+'Obj'+pos.OBJID+'P'+pos.POS+'_slitlineI_'+linetag+'_bsw'+str(slitpixsize)+'ars.fits'
	filein_fnu=dir_obj+'Obj'+pos.OBJID+'P'+pos.POS+'_slitfnu'+'_bsw'+str(slitpixsize)+'ars.fits'

	z=pos.outer.outer.ZSTAR
	bandratiopairs=[['r','u'],['r','g'],['r','i'],['r','z']]
	extcolorpairs=[['u','r'],['r','z']]
	noisethreshold=0.02 # nanomaggy
	seterr(all='ignore')

	# read in
	if os.path.isfile(filein_lI):
		bintable=Table.read(filein_lI,format='fits')
		bintable=join(bintable,Table.read(filein_fnu,format='fits'))
	else: 
		print "skip incorporating I_lOIII5008 in to SlitAll Table as input file does not exist"
		bintable=Table.read(filein_fnu,format='fits')

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




def poswSlitFluxDensity(pos,bands=['u','g','r','i','z'],imagepixsize=0.396,slitpixsize=0.6,slitlength=20.,slitwidth=1.0,toplot=True):
	"""
	PURPOSE: obtain flux density along magellan slit position given "pos" object

	DESCRIPTION: 
		Flux density values extracted from averaged spline interpolated image values at positions around the pixel location along the slit crossection with width=slitwidth. 

	INPUT: 
		pos: position object
		imagepixsize: pix size of image [arcsec]
		slitpixsize:  pix size of slit  [arcsec]
		slitlength:   length of slit    [arcsec]
		slitwidth:    width of slit     [arcsec]


	WRITE OUTPUT:  
		e.g., Obj2100P1_slitfnu_bsw0.6ars.fits

		content: bintable
		binindex  ycoord    Fnu_g_nmaggy      Fnu_r_nmaggy    Fnu_i_nmaggy  
		          arcsec                                                    
		 int64   float64      float64           float64         float64     
		-------- ------- ------------------ --------------- ----------------

	PLOT OUTPUT: 
		e.g., Obj2100P1_slitfnu_bsw0.6ars.pdf
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
		filein=dir_obj+'stamp-'+band+'.fits'
		hdulist=fits.open(filein)
		image=hdulist[0].data
		header=hdulist[0].header
		# image spec
		nx,ny=image.shape
		x0,y0=header['CRPIX1'],header['CRPIX2']
		# spline 
		f=interpolate.interp2d(arange(nx), arange(ny), image, kind='linear')
		# relative PA (alpha)
		imagePA=getimagePA_from_header(hdulist[0].header)
		alpha=pos.PA-imagePA

		for i in range(nbin):
			r=bintable['ycoord'][i]
			# # (x,y) is a vector from image center to the bin pixel centroid in units of image pixel
			xl, yl, xr, yr=getlineedges(r,alpha,slitwidth,imagepixsize,x0,y0)

			# simulate slit with finite width by sampling linearly interpolated image with a few points in between (xl,yl) and (xr,yr)
			xs=linspace(xl,xr,num=10.)
			ys=linspace(yl,yr,num=10.)
			
			# record slit pixel value as the average of the interpolated image values around the bin pixel
			bintable['Fnu_'+band+'_nmaggy'][i]=average([f(xs[j],ys[j])[0]  for j in range(len(xs))])

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



def poswSlitLineIntensity(pos,linetag='lOIII5008',imagepixsize=0.396,slitpixsize=0.6,slitlength=20.,slitwidth=1.0,toplot=True):
	"""
	PURPOSE: obtain line intensity [erg/cm2/s/arcsec2] along magellan slit position given "pos" object

	DESCRIPTION: 
		Line intensity values extracted from averaged spline interpolated image values at positions around the pixel location along the slit crossection with width=slitwidth. 

	INPUT: 
		pos: position object
		linetag='lOIII5008'
		imagepixsize: pix size of image [arcsec]
		slitpixsize:  pix size of slit  [arcsec]
		slitlength:   length of slit    [arcsec]
		slitwidth:    width of slit     [arcsec]
		toplot=True


	WRITE OUTPUT:  
		e.g., Obj2100P1_slitlineI_lOIII5008_bsw0.6ars.fits

		content: bintable
		binindex  ycoord    I_lOIII5008      
		          arcsec  erg/cm2/s/arcsec2                                                
		 int64   float64      float64       
		-------- ------- ------------------ 

	PLOT OUTPUT: 
		e.g., Obj2100P1_slitlineI_lOIII5008_bsw0.6ars.pdf
	"""

	# setup
	dir_obj=dir_data+'M'+pos.OBJID+'/'
	fileout=dir_obj+'Obj'+pos.OBJID+'P'+pos.POS+'_slitlineI_'+linetag+'_bsw'+str(slitpixsize)+'ars'

	# read in
	filein=dir_obj+'stamp-'+linetag+'_I'+'.fits'
	if not os.path.isfile(filein): print 'skip poswSlitLineIntensity as input file stamp-'+linetag+'_I'+'.fits does not exist'
	else:
		hdulist=fits.open(filein)
		image=hdulist[0].data
		header=hdulist[0].header

		# initialize bintable, containing columns ("binindex" [pix]  "ycoord" [arcsec])
		nbin=int(floor(slitlength/slitpixsize / 2.) * 2 +1) # number of bins rounded to closest odd number
		binindexes=arange(nbin)-(nbin-1)/2
		binycoord=binindexes*slitpixsize
		bintable=Table(data=[binindexes,binycoord],names=['binindex','ycoord'])
		bintable['ycoord'].unit=u.arcsec

		# writing line intensity
		bintable['I_'+linetag]=0.
		bintable['I_'+linetag].unit=header['BUNIT']

		# image spec
		nx,ny=image.shape
		x0,y0=header['CRPIX1'],header['CRPIX2']
		# spline interpolate image
		f=interpolate.interp2d(arange(nx), arange(ny), image, kind='linear')
		# relative PA (alpha)
		imagePA=getimagePA_from_header(hdulist[0].header)
		alpha=pos.PA-imagePA

		for i in range(nbin):
			r=bintable['ycoord'][i]
			xl, yl, xr, yr=getlineedges(r,alpha,slitwidth,imagepixsize,x0,y0)

			# simulate slit with finite width by sampling linearly interpolated image with a few points in between (xl,yl) and (xr,yr)
			xs=linspace(xl,xr,num=10.)
			ys=linspace(yl,yr,num=10.)
			# plot(xs,ys,color='red')

			# record slit pixel value as the average of the interpolated image values around the bin pixel
			bintable['I_'+linetag][i]=average([f(xs[j],ys[j])[0]  for j in range(len(xs))])

		# output
		bintable.write(fileout+'.fits',format='fits',overwrite=True)
		bintable.write(fileout+'.txt',format='ascii.fixed_width',delimiter='')

		if toplot:
			clf() # plot line intensity profile
			plot(bintable['ycoord'],bintable['I_'+linetag],label='I_'+linetag)
			axhline(0,color='black')
			legend()
			xlabel('ycoord ['+bintable['ycoord'].unit.to_string()+']')
			ylabel('I ['+bintable['I_'+linetag].unit.to_string()+']')
			savefig(fileout+'.pdf')

			clf() # plot slit position on top of SDSS image
			imshow(image,origin='lower',cmap='gist_gray',vmin=0.,vmax=image.max()*0.5,interpolation='nearest')
			rs=bintable['ycoord']

			xls, yls, xrs, yrs=getlineedges(rs,alpha,slitwidth,imagepixsize,x0,y0)

			plot(xls,yls,ls='-',color='white',lw=2)
			plot(xrs,yrs,ls='-',color='white',lw=2)
			xlim(0,nx-1)
			ylim(0,ny-1)
			savefig(fileout+'_slitpos.pdf')



def getlineedges(rs,alpha,slitwidth,imagepixsize,x0,y0):
	"""
	Return line edges positions (xls, yls) for the left line, and (xrs, yrs) for the right line, given:
		rs: array of radial distance(s) from origin (x0, y0) [arcsec]
		alpha: angle of the slit relative to image, pos.PA-imagePA
		slitwidth: [arcsec]
		imagepixsize: [arcsec]
		x0,y0: origin or reference poitn on the image [pix]

	It's an advance version of marking the slit centroid line
			# x=-r*sin(radians(alpha))/imagepixsize+x0 
			# y= r*cos(radians(alpha))/imagepixsize+y0

	"""
	# print alpha
	xls=(-rs*sin(radians(alpha))+0.5*slitwidth*cos(radians(alpha)))/imagepixsize+x0
	yls=( rs*cos(radians(alpha))+0.5*slitwidth*sin(radians(alpha)))/imagepixsize+y0
	xrs=(-rs*sin(radians(alpha))-0.5*slitwidth*cos(radians(alpha)))/imagepixsize+x0
	yrs=( rs*cos(radians(alpha))-0.5*slitwidth*sin(radians(alpha)))/imagepixsize+y0
	return xls, yls, xrs, yrs


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
		pa=180.+degrees(arcsin(sintheta))

	return pa



