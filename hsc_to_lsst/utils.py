from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
import warnings
import numpy as np
from scipy.interpolate import UnivariateSpline


def photutils_background_iterative(
        data,
        nsigma_detection=3,
        sigma_clip=3,
        clip_iters=10,
        npixels_detection=5,
        init_median=0,
        init_rms=None,
        mask_size=2,
        iters=3
):
    if init_rms is None:
        init_rms = sigma_clipped_stats(data, sigma=sigma_clip, maxiters=clip_iters)[-1]
    median, std = init_median, init_rms
    mask = np.zeros_like(data, dtype=bool)
    for i in range(iters):
        threshold = detect_threshold(data, nsigma=nsigma_detection, background=median, error=std)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            segment_img = detect_sources(data, threshold, npixels=npixels_detection)
        footprint = circular_footprint(radius=mask_size)
        if segment_img is None:
            break
        mask = segment_img.make_source_mask(footprint=footprint)
        if np.all(mask):
            break
        mean, median, std = sigma_clipped_stats(data, sigma=sigma_clip, mask=mask, maxiters=clip_iters)
    return median, std, mask


def get_fwhm(psf, pix_scale):
    # Define the center of the PSF (assuming it's approximately centered)
    center_x = (psf.shape[1] - 1) // 2
    center_y = (psf.shape[0] - 1) // 2

    # Extract the row and column passing through the center
    central_row = psf[center_y, :]
    central_col = psf[:, center_x]

    half_max = np.max(psf) / 2

    # Interpolating the central row and column to find FWHM
    spline_row = UnivariateSpline(np.arange(len(central_row)), central_row - half_max, s=0)
    spline_col = UnivariateSpline(np.arange(len(central_col)), central_col - half_max, s=0)

    # Find the points where the profile crosses half maximum
    fwhm_row = np.abs(spline_row.roots()[0] - spline_row.roots()[-1])
    fwhm_col = np.abs(spline_col.roots()[0] - spline_col.roots()[-1])

    # Average FWHM (for a circular PSF)
    fwhm = (fwhm_row + fwhm_col) / 2

    return fwhm * pix_scale
