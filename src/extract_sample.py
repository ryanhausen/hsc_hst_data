import os
import random
from functools import partial
from itertools import takewhile
from typing import Tuple

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from tqdm import tqdm
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
    ys, xs = slice(y-size//2, y+size//2), slice(x-size//2, x+size//2)
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


    try:
        hst_sample = Cutout2D(
            data=hst.data,
            position=skycoord,
            size=hst_size,
            wcs=WCS(hst.header),
        )
    except NoOverlapError:
        print(f"x coord: {x}, y coord: {y}")


    return (
        fits.PrimaryHDU(data=hsc_sample.data, header=hsc_sample.wcs.to_header()),
        fits.PrimaryHDU(data=hst_sample.data, header=hst_sample.wcs.to_header()),
    )

def random_sample_generator(y_bnds, x_bnds):
    while True:
        yield random.randint(*y_bnds), random.randint(*x_bnds)

def main():

    mask = fits.getdata("./data/hst_to_hsc_footprint.fits")
    hsc = fits.open("./data/cutout-HSC-I-9813-pdr2_dud-210317-161628.fits")[1]
    hst = fits.open("./data/hlsp_candels_hst_acs_cos-tot_f814w_v1.0_drz.fits")[0]

    hsc_sample_size = 100
    hst_sample_size = hsc_sample_size * HST_HSC_RATIO

    # ==========================================================================
    # Add pixel locations here!
    # ========= =================================================================
    validate_idx_f = partial(validate_sample, mask, hsc_sample_size)

    num_samples = 10 # Number of samples to try to generate
    edge_scaler = hsc_sample_size/2    # Handle edges
    y_max,x_max =[dim-edge_scaler for dim in hsc.shape] # Subtract 1/2 dimension to avoid edges
    y_bounds, x_bounds = (edge_scaler, y_max), (edge_scaler, x_max)

    hsc_sample_locations = takewhile(
        lambda idx_coords: idx_coords[0] < num_samples,
        enumerate(
            filter(
                lambda yx: validate_idx_f(*yx),
                random_sample_generator(y_bounds, x_bounds)
            )
        )
    )

    extract_function_f = partial(
        extract_sample,
        hsc,
        hst,
        hsc_sample_size,
        hst_sample_size
    )

    data_dirs = [
        "./data/samples/hsc",
        "./data/samples/hst"
    ]

    for dd in data_dirs:
        if not os.path.exists(dd):
            os.makedirs(dd)

    for idx, yx in tqdm(hsc_sample_locations, total=num_samples):
        hsc_sample, hst_sample = extract_function_f(*yx)
        # UNCOMMENT BELOW TO GENERATE SAMPLES
        hsc_sample.writeto(f"./data/samples/hsc/{idx}.fits", overwrite=True)
        hst_sample.writeto(f"./data/samples/hst/{idx}.fits", overwrite=True)




if __name__=="__main__":
    main()
