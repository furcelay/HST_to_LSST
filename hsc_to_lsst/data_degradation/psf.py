import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from scipy.interpolate import UnivariateSpline


def psf_kernel_from_fhwm_diff(original_fwhm, new_fwhm, original_pix_scale):
    trans_fwhm = np.sqrt(new_fwhm**2 - original_fwhm**2)
    kernel_sigma = trans_fwhm * gaussian_fwhm_to_sigma / original_pix_scale
    kernel = Gaussian2DKernel(x_stddev=kernel_sigma)
    return kernel


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
