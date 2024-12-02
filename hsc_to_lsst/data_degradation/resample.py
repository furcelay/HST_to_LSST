from astropy.wcs import WCS
from reproject import reproject_adaptive
import numpy as np


deg2arcsec = 3600


def resample_image(img, original_pix_scale, lsst_pix_scale=0.2, drop_edge=5, out_size=None):
    img_wcs = WCS(naxis=2)
    img_wcs.wcs.cd = np.array([[-1, 0],
                               [ 0, 1]]) * original_pix_scale / deg2arcsec
    img_crpix = np.flip(np.array(img.shape)) / 2
    img_wcs.wcs.crpix = img_crpix

    img_field = np.array(img.shape) * original_pix_scale
    img_center = img_wcs.pixel_to_world(*((np.array(img.shape) - 1) / 2))

    if out_size is None:
        lsst_size = img_field // lsst_pix_scale
    else:
        lsst_size = out_size
    lsst_crpix = np.flip(lsst_size) / 2

    lsst_wcs = WCS(naxis=2)
    lsst_wcs.wcs.crpix = lsst_crpix
    lsst_wcs.wcs.cd = np.array([[-1, 0],
                                [ 0, 1]]) * lsst_pix_scale / deg2arcsec
    lsst_wcs.wcs.ctype = img_wcs.wcs.ctype

    img_scaled, _ = reproject_adaptive((img, img_wcs), lsst_wcs,
                                       shape_out=lsst_size.astype(int), conserve_flux=True,
                                       bad_fill_value=0)
    # drop edge that has bad pixels
    if drop_edge > 0:
        img_scaled = img_scaled[drop_edge:-drop_edge,
                                drop_edge:-drop_edge]
    return img_scaled