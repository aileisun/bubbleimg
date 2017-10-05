# hscurl.py
# ALS 2017/06/28


def get_hsc_cutout_url(ra, dec, band='i', rerun='', tract='', imgtype='coadd', sw='5asec', sh='5asec', mask='on', variance='on'):
	"""
	see hsc query manual
	https://hscdata.mtk.nao.ac.jp/das_quarry/manual.html 
	"""

	# public
	# https://hsc-release.mtk.nao.ac.jp/das_quarry/cgi-bin/quarryImage?ra=-24&dec=0&sw=2asec&sh=2asec&type=coadd&image=on&filter=HSC-G&tract=&rerun=

	# old
	# url = 'https://hscdata.mtk.nao.ac.jp:4443/das_quarry/cgi-bin/quarryImage?ra={0}&dec={1}&sw={2}&sh={3}&type={4}&image=on&mask={5}&variance={6}&filter=HSC-{7}&tract={8}&rerun={9}'.format(ra, dec, sw, sh, imgtype, mask, variance, band.capitalize(), tract, rerun)

	url = 'https://hscdata.mtk.nao.ac.jp/das_quarry/dr1/cgi-bin/quarryImage?ra={0}&dec={1}&sw={2}&sh={3}&type={4}&image=on&mask={5}&variance={6}&filter=HSC-{7}&tract={8}&rerun={9}'.format(ra, dec, sw, sh, imgtype, mask, variance, band.capitalize(), tract, rerun)

	return url


def get_hsc_psf_url(ra, dec, band='i', rerun='', tract='', patch='', imgtype='coadd'):
	"""
	see hsc query manual
	https://hscdata.mtk.nao.ac.jp/psf/4/manual.html#Bulk_mode
	"""

	url = 'https://hscdata.mtk.nao.ac.jp/psf/4/cgi/getpsf?ra={ra}&dec={dec}&filter={band}&rerun={rerun}&tract={tract}&patch={patch}&type={imgtype}'.format(ra=ra, dec=dec, band=band, rerun=rerun, tract=tract, patch=patch, imgtype=imgtype)
	return url

