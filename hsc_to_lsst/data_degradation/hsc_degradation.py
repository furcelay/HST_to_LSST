from astropy.nddata import Cutout2D
from lenstronomy.SimulationAPI.ObservationConfig.LSST import LSST
from hsc_to_lsst.data_degradation.psf import degrade_psf
from hsc_to_lsst.data_degradation.resample import resample_image
from hsc_to_lsst.data_degradation.noise import add_noise


def hsc_to_lsst(
        hsc_img,
        lsst_band,
        exp_time=30.0,
        lsst_zero_point=31.0,
        hsc_fwhm=0.6,
        hsc_psf=None,
        hsc_pix_scale=0.168,
        lsst_fwhm=None,
        lsst_psf=None,
        background_noise=None,
        background_median=0,
        psf_transform=True,
        add_poisson_noise=True,
        add_background_noise=True,
        use_nise_diff=True,
        to_adu=False,
        out_size=64
):
    lsst_band_props = LSST(band=lsst_band).kwargs_single_band()
    if lsst_fwhm is None:
        lsst_fwhm = lsst_band_props['seeing']
    lsst_pix_scale = lsst_band_props['pixel_scale']
    # PSF
    if psf_transform:
        hsc_img_conv = degrade_psf(
                hsc_img,
                original_fwhm=hsc_fwhm,
                target_fwhm=lsst_fwhm,
                original_psf=hsc_psf,
                target_psf=lsst_psf,
                original_pix_scale=hsc_pix_scale,
                target_pix_scale=lsst_pix_scale,
                max_iters=3,
                thresh=0.01
            )
    else:
        hsc_img_conv = hsc_img

    # resampling
    hsc_img_conv_scaled = resample_image(hsc_img_conv, hsc_pix_scale, lsst_pix_scale, drop_edge=0)
    # cut to out size, centered
    cy = (hsc_img_conv_scaled.shape[0] - 1) / 2
    cx = (hsc_img_conv_scaled.shape[1] - 1) / 2
    hsc_img_conv_scaled = Cutout2D(hsc_img_conv_scaled, (cy, cx), (out_size, out_size)).data
    # noise
    if add_poisson_noise or add_background_noise:
        hsc_img_conv_scaled_noise = add_noise(hsc_img_conv_scaled,
                                              lsst_band_props,
                                              exp_time,
                                              lsst_zero_point,
                                              background_noise=background_noise,
                                              background_median=background_median,
                                              add_poisson_noise=add_poisson_noise,
                                              add_background_noise=add_background_noise,
                                              use_nise_diff=use_nise_diff
                                              )
    else:
        hsc_img_conv_scaled_noise = hsc_img_conv_scaled
    if to_adu:
        hsc_img_conv_scaled_noise /= lsst_band_props['ccd_gain']

    return hsc_img_conv_scaled_noise
