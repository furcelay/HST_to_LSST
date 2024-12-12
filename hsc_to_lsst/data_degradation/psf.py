import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from hsc_to_lsst.utils import get_fwhm
from astropy.convolution import convolve
import warnings


def psf_kernel_from_fhwm(fwhm, pix_scale):
    kernel_sigma = fwhm * gaussian_fwhm_to_sigma / pix_scale
    kernel = Gaussian2DKernel(x_stddev=kernel_sigma)
    return kernel


def psf_kernel_from_fhwm_diff(original_fwhm, target_fwhm, original_pix_scale):
    trans_fwhm = np.sqrt(target_fwhm**2 - original_fwhm**2)
    return psf_kernel_from_fhwm(trans_fwhm, original_pix_scale)


def iterative_psf_transform_kernel(original_psf, target_fwhm,
                                   original_pix_scale,
                                   max_iters=3, thresh=0.01):
    original_fwhm = get_fwhm(original_psf, original_pix_scale)
    correction_factor = 0
    trans_kernel = psf_kernel_from_fhwm_diff(original_fwhm + correction_factor, target_fwhm, original_pix_scale)
    for i in range(max_iters):
        transformed_psf = convolve(original_psf, trans_kernel)
        transformed_fwhm = get_fwhm(transformed_psf, original_pix_scale)
        if abs(transformed_fwhm - target_fwhm) < thresh:
            break
        correction_factor += transformed_fwhm - target_fwhm
        trans_kernel = psf_kernel_from_fhwm_diff(original_fwhm + correction_factor, target_fwhm, original_pix_scale)
    return trans_kernel


def degrade_psf(
        image,
        original_fwhm=None,
        target_fwhm=None,
        original_psf=None,
        target_psf=None,
        original_pix_scale=0.168,
        target_pix_scale=0.20,
        max_iters=3,
        thresh=0.01
):
    if target_psf is not None:
        target_fwhm = get_fwhm(target_psf, target_pix_scale)
    if original_psf is not None:
        original_fwhm = get_fwhm(original_psf, original_pix_scale)
        if original_fwhm > target_fwhm:
            warnings.warn("Warning: original FWHM is larger than target")
            return image
        trans_kernel = iterative_psf_transform_kernel(original_psf, target_fwhm, original_pix_scale, max_iters, thresh)
    else:
        if original_fwhm > target_fwhm:
            warnings.warn("Warning: original FWHM is larger than target")
            return image
        trans_kernel = psf_kernel_from_fhwm_diff(original_fwhm, target_fwhm, original_pix_scale)
    return convolve(image, trans_kernel)
