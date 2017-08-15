# plotpsf_tools.py

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def _normalize_trim_psf(psf0, psf1, psfmt, w=21):
	psf0 = psf0/psf0.sum()
	psf1 = psf1/psf1.sum()
	psfmt = psfmt/psfmt.sum()

	pmin = min([psf0.min(), psf1.min(), psfmt.min()])
	pmax = max([psf0.max(), psf1.max(), psfmt.max()])

	psf0 = psf0/pmax
	psf1 = psf1/pmax
	psfmt = psfmt/pmax

	psf0 = _cut_psf_center(psf0, w=w)
	psf1 = _cut_psf_center(psf1, w=w)
	psfmt = _cut_psf_center(psfmt, w=w)

	return psf0, psf1, psfmt


def _cut_psf_center(psf, w=21):
	nx, ny = psf.shape

	x0 = (nx-w)/2
	x1 = nx - x0

	y0 = (ny-w)/2
	y1 = ny - y0

	psf_cut = psf[x0:x1, y0:y1]

	return psf_cut


def _plot_psfmatch_kernel(band, bandto, psf0, psf1, psfmt, fn):
	vmin = 0
	vmax = 1

	f = plt.figure(figsize=(16, 5))
	f.clf()

	gs = gridspec.GridSpec(2, 4, height_ratios=[1, 0.05, ], hspace=0.0, wspace=0.2, left=0.05, right=0.95)

	ax1 = plt.subplot(gs[0, 0])
	ax2 = plt.subplot(gs[0, 1])
	ax3 = plt.subplot(gs[0, 2])
	ax4 = plt.subplot(gs[0, 3])

	# cx1 = plt.subplot(gs[1, 0])
	# cx2 = plt.subplot(gs[1, 1])
	cx3 = plt.subplot(gs[1, 2])
	cx4 = plt.subplot(gs[1, 3])

	im = ax1.imshow(psf0, vmin=vmin, vmax=vmax)
	ax1.set_title('(a) {}-band'.format(band), fontsize=18)
	# f.colorbar(im, cax=cx1, orientation='horizontal')

	im = ax2.imshow(psf1, vmin=vmin, vmax=vmax)
	ax2.set_title('(b) {}-band'.format(bandto), fontsize=18)
	# f.colorbar(im, cax=cx2, orientation='horizontal')

	im = ax3.imshow(psfmt, vmin=vmin, vmax=vmax)
	ax3.set_title('(c) {}-band matched to {}-band'.format(band, bandto), fontsize=18)
	f.colorbar(im, cax=cx3, orientation='horizontal')

	im = ax4.imshow(psfmt-psf1)
	ax4.set_title('(d) difference between (c) and (b)', fontsize=18)
	cb = f.colorbar(im, cax=cx4, orientation='horizontal')

	cb.formatter.set_powerlimits((0, 0))
	cb.update_ticks()

	# remove the x and y ticks
	for ax in (ax1, ax2, ax3, ax4):
		ax.set_xticks([])
		ax.set_yticks([])

	f.savefig(fn)
	plt.close(f)
