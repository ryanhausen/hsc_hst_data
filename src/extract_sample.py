import os
from functools import partial
from typing import Tuple

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D

# HST arcsec/pixel = 8.333333E-6
# HSC arcsec/pixel = 4.66666666666396E-05

HST_HSC_RATIO = 4.66666666666396E-05 / 8.333333E-6
#WCS_HSC = WCS(fits.open("../data/cutout-HSC-I-9813-pdr2_dud-210317-161628.fits")[1].header)
#WCS_HST =  WCS("../data/hlsp_candels_hst_acs_cos-tot_f814w_v1.0_drz.fits")


def validate_sample(mask:np.ndarray, size:int, y:int, x:int) -> bool:
    """Validates a coordinate(y,x) and size for HSC/HST compatability.

    args:
        mask (np.ndarray): A boolean array that where True indicates that a
                           pixel exists in both HSC/HST and False indicates
                           it exists in the HSC only.
        y (int): y coordinate in image space representing sample center.
        x (int): x coordinate in image space representing sample center.

    returns
        True if the size^2 sample at (y, x) exists in both HST and HSC and False
        if it only exists in one of the images.
    """
    ys, xs = slice(y-size, y+size), slice(x-size, x+size)
    return mask[ys, xs].all()

def extract_sample(
    hsc:fits.PrimaryHDU,
    hst:fits.PrimaryHDU,
    hsc_size:int,
    hst_size:int,
    y:int,
    x:int,
) -> Tuple[fits.PrimaryHDU, fits.PrimaryHDU]:
    """Extracts samples from the HSC/HST images that correpsond to the same area of the sky.

    args:
        hsc (fits.PrimaryHDU): The HSC HDU to extract the sample from.
        hst  (fits.PrimaryHDU): The HST HDU to extract the sample from.
        hsc_size (int): The size of the HSC cutout where width=height=size
        hst_size (int): The size of the HST cutout where width=height=size
        y (int): y coordinate in image space representing sample center.
        x (int): x coordinate in image space representing sample center.

    returns
        A tuple where the first element is the HSC sample and the second element
        is the HST sample.
    """
    skycoord = WCS(hsc.header).pixel_to_world(x, y)

    hsc_sample = Cutout2D(
        data=hsc.data,
        position=skycoord,
        size=hsc_size,
        wcs=WCS(hsc.header)
    )

    hst_sample = Cutout2D(
        data=hst.data,
        position=skycoord,
        size=hst_size,
        wcs=WCS(hst.header),
    )

    return (
        fits.PrimaryHDU(data=hsc_sample.data, header=hsc_sample.wcs.to_header()),
        fits.PrimaryHDU(data=hst_sample.data, header=hst_sample.wcs.to_header()),
    )

def extract_sample_old(y, x):
    mask = fits.getdata("../data/hst_to_hsc_footprint.fits")

    hsc_sample_y, hsc_sample_x = 10229, 4383
    hsc_sky_coord = WCS_HSC.pixel_to_world(hsc_sample_x, hsc_sample_y)
    hsc_size = 100
    hsc = fits.getdata("../data/cutout-HSC-I-9813-pdr2_dud-210317-161628.fits")
    hsc_sample = Cutout2D(hsc, hsc_sky_coord, hsc_size, wcs=WCS_HSC)

    fits.PrimaryHDU(data=hsc_sample.data, header=hsc_sample.wcs.to_header()).writeto("sample_hsc.fits", overwrite=True)


    hst_size = hsc_size * HST_HSC_RATIO
    hst = fits.getdata("../data/hlsp_candels_hst_acs_cos-tot_f814w_v1.0_drz.fits")
    hst_sample = Cutout2D(hst, hsc_sky_coord, hst_size, wcs=WCS_HST)

    fits.PrimaryHDU(data=hst_sample.data, header=hst_sample.wcs.to_header()).writeto("sample_hst.fits", overwrite=True)

def main():

    mask = fits.getdata("./data/hst_to_hsc_footprint.fits")
    hsc = fits.open("./data/cutout-HSC-I-9813-pdr2_dud-210317-161628.fits")[1]
    hst = fits.open("./data/hlsp_candels_hst_acs_cos-tot_f814w_v1.0_drz.fits")[0]

    hsc_sample_size = 100
    hst_sample_size = hsc_sample_size * HST_HSC_RATIO

    # ==========================================================================
    # Add pixel locations here!
    # ==========================================================================
    hsc_sample_locations = [
        (10229, 4383)
    ]

    validate_idx_f = partial(validate_sample, mask, hsc_sample_size)
    extract_function_f = partial(
        extract_sample,
        hsc,
        hst,
        hsc_sample_size,
        hst_sample_size
    )

    data_dirs = [
        "./data/samples/hsc"
        "./data/samples/hst"
    ]
    for dd in data_dirs:
        if not os.path.exists(dd):
            os.path.makedirs(dd)

    for idx, yx in enumerate(hsc_sample_locations):
        if validate_idx_f(*yx):
            hsc_sample, hst_sample = extract_function_f(*yx)

            hsc_sample.writeto(f"./data/samples/hsc/{idx}.fits", overwrite=True)
            hst_sample.writeto(f"./data/samples/hst/{idx}.fits", overwrite=True)




if __name__=="__main__":
    main()
