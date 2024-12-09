from hsc_to_lsst.utils import photutils_background_iterative
from lenstronomy.Util import data_util
import numpy as np
# import warnings


def add_noise(img,
              lsst_band_props,
              exp_time=30.0,
              zero_point=31.0,
              background_noise=None,
              background_median=0,
              add_poisson_noise=True,
              add_background_noise=True,
              use_nise_diff=True,
              ):
    img_median, img_std, _ = photutils_background_iterative(img)
    img = img - img_median

    num_exposures = exp_time / 15

    img_positive = np.where(img > 0, img, 0)

    # image with poisson noise in e-/s
    if add_poisson_noise:
        img = np.random.poisson(lam=img_positive * exp_time) / exp_time

    if add_background_noise:
        if background_noise is None:
            sky_brightness_cps = data_util.magnitude2cps(
                    lsst_band_props['sky_brightness'],
                    magnitude_zero_point=zero_point,
                )

            bkg_noise = data_util.bkg_noise(
                    lsst_band_props['read_noise'],
                    15,
                    sky_brightness_cps,
                    lsst_band_props['pixel_scale'],
                    num_exposures,
                )
        else:
            bkg_noise = background_noise
        if bkg_noise > img_std:
            if use_nise_diff:
                noise = np.random.normal(scale=np.sqrt(bkg_noise**2 - img_std**2),
                                         size=img.shape)
                img = img + noise
            else:
                noise = np.random.normal(scale=bkg_noise, size=img.shape)
                img = img + noise
        else:
            # warnings.warn("Warning: original noise is larger than target")
            raise ValueError("Background noise is larger than target")

    img_median = photutils_background_iterative(img)[0]
    img += background_median - img_median

    return img
