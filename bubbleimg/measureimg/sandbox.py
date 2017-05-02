# sandbox.py

def add_gaussnoise(img, noise_sigma):
    """ add iid gaussian noise on image"""
    from skimage.util import random_noise
    if noise_sigma==0:
        return img 
    else: 
        return random_noise(img, mode='gaussian',var=noise_sigma**2,clip=False)



nx,ny=64,64
noise_sigma=1.
img=np.zeros([nx,ny])

img=add_gaussnoise(img, noise_sigma)

# img[31:33,31:33]=10.

threshold=1.
areallimit=3.
# contours=find_isocontours(imgn,threshold)
contours=find_realisocontours(img, threshold, areallimit)

print len(contours)
clf()
imshow(img)
for contour in contours:
    plot(contour[:, 1], contour[:, 0],color='black')


def sandbox():
    poly=maincontour
    plt.close('all')
    plt.plot(poly[0:21,0],poly[0:21,1])
    plt.plot(poly[20:41,0],poly[20:41,1])
    plt.plot(poly[40:61,0],poly[40:61,1])
    plt.plot(poly[60:81,0],poly[60:81,1])
    plt.plot(poly[80:101,0],poly[80:101,1])
    plt.plot(poly[100:121,0],poly[100:121,1])
    plt.show(block=False)


filename='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/algorithm/bubblepy/testing/SDSSJ1000+1242/measureimg_iso_stamp-lOIII5008_I_norm_denoised_contours.pkl'
with open(filename, 'rb') as handle:
    b=pickle.load(handle)

contrast=0.1

filename_gal='stamp-conti-onOIIIscale_I_norm.fits'
isocut_rest=1./contrast*3.e-15*u.Unit('erg / (arcsec2 cm2 s)')
isoareallimit=0.

contours_galshade, isothreshold = dir_makeContours(dir_obj, filename=filename_gal, isocut_rest=isocut_rest, isoareallimit=isoareallimit, update=False)

filename='stamp-lOIII5008_I_norm_denoised.fits'
img=fits.getdata(dir_obj+filename)
isocut_rest=3.e-15*u.Unit('erg / (arcsec2 cm2 s)')
isoareallimit=10.
contours, isothreshold = dir_makeContours(dir_obj, filename=filename, isocut_rest=isocut_rest, isoareallimit=isoareallimit, update=False)


close()
fig,ax=plottools.plot_img(img,vmin=None,vmax=None)
plottools.overplot_contours(ax, contours_galshade, color='black')
plottools.overplot_contours(ax, contours, color='white')




contours

contours_galshade

contours2=[np.array([[0,28],[63,28],[63,32],[0,32],[0,28]])]
# contours2=[np.array([[2,28],[2,32],[63,32],[63,28],[2,28]])]

sgM=sgMultipolygon_from_contours(contours)
sgM_galshade=sgMultipolygon_from_contours(contours_galshade)

sgMdiff=sgM.difference(sgM_galshade)


contoursdiff=diffcontours(contours, contours_galshade)

plottools.overplot_contours(ax, contoursdiff, color='red')


contoursdiff2=diffcontours(contours, contours2)
plottools.overplot_contours(ax, contours2, color='green')
plottools.overplot_contours(ax, contoursdiff2, color='red')


with open(a, 'rb') as handle:
    cdict_linblobmasked=pickle.load(handle)