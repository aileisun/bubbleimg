# standards.py
"""
setting standards for the pipeline
"""


def get_img_xycenter(img):
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
    return get_img_xycenter_fromnxny(nx, ny)


def get_img_xycenter_fromnxny(nx, ny):
    """ see get_img_xycenter """
    if (nx % 2 == 0) and (ny % 2 == 0):
        xc, yc = nx/2, ny/2
    elif (nx % 2 == 1) and (ny % 2 == 1):
        xc, yc = (nx-1)/2, (ny-1)/2
    else:
        raise ValueError("img shape is not accepted")

    return xc, yc
