# standards.py
"""
setting standards for the pipeline
"""


def get_img_xycenter(img, center_mode='n/2'):
    """
    use the convention of alignshift to define img center.
    If nx, ny are even, the center is on [nx/2., ny/2.]
    Otherwise if odd, the center is on [(nx-1)/2.,(ny-1)/2.]

    Params
    --------
    img: np 2d array

    Return
    --------
    xc: int
    yc: int
    """
    nx, ny = img.shape

    xc, yc = get_img_xycenter_fromnxny(nx, ny, center_mode=center_mode)

    return xc, yc


def get_img_xycenter_fromnxny(nx, ny, center_mode='n/2'):
    """ see get_img_xycenter """

    xc = get_center_from_n(nx, center_mode=center_mode)
    yc = get_center_from_n(ny, center_mode=center_mode)

    return xc, yc


def get_center_from_n(n, center_mode='n/2'):
    if (n % 2 == 0):
        if center_mode == 'n/2':
            c = n/2
        elif center_mode == 'n/2-1':
            c = n/2-1
        else: 
            raise ValueError("[standards] param center_mode not recognized")
    elif (n % 2 == 1):
        c = (n-1)/2

    return c

