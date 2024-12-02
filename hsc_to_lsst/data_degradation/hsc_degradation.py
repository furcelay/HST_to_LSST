def hsc_to_lsst(hsc_img,
                hsc_hdr,
                lsst_band,
                num_visits=1,
                exp_time=30.0,
                lsst_zero_point=31.0,
                # max_zero_point_diff=3,
                hsc_fwhm=0.6,
                lsst_fwhm=None,
                background_noise=None,
                background_median=0,
                # change_zero_point=True,
                psf_transform=True,
                add_poisson_noise=True,
                add_background_noise=True,
                use_nise_diff=True,
                to_adu=False,
                out_size=64):
    hsc_wcs = WCS(hsc_hdr)
    try:
        hsc_pix_scale = hsc_wcs.wcs.cd[1, 1] * deg2arcsec
    except AttributeError:
        hsc_pix_scale = hsc_wcs.wcs.pc[1, 1] * deg2arcsec
    lsst_band_props = LSST(band=lsst_band).kwargs_single_band()
    if lsst_fwhm is None:
        lsst_fwhm = lsst_band_props['seeing']
    lsst_pix_scale = lsst_band_props['pixel_scale']
    # Zero Point
    # if change_zero_point:
    #     hsc_img = zero_point_change(hsc_img, 27.0, lsst_zero_point, max_zero_point_diff)
    # PSF
    if psf_transform:
        if lsst_fwhm > hsc_fwhm:
            lsst_kernel = psf_kernel_from_fhwm_diff(hsc_fwhm, lsst_fwhm, hsc_pix_scale)
            hsc_img_conv = convolve(hsc_img, lsst_kernel, boundary='extend')
        else:
            warnings.warn("Warning: original FWHM is larger than target")
            hsc_img_conv = hsc_img
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
                                            num_visits,
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