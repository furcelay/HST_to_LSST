from hsc_to_lsst.hsc_query import query_hsc
from hsc_to_lsst.data_degradation.zero_point import zero_point_change
from hsc_to_lsst.data_degradation.hsc_degradation import hsc_to_lsst
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import warnings


HSC_FWHM = {
    'g': 0.71,
    'r': 0.80,
    'i': 0.53,
    'z': 0.56,
    'y': 0.53
}


def query_and_degrade(
        ra,
        dec,
        size,
        field,
        username,
        password,
        dp0_sampler,
        num_visits=1,
        zp_rms_frac_thresh=0.1
):
    try:
        hsc_data = query_hsc(username, password, ra, dec, size, field)
    except OSError:
        return False, None, None

    pix_scale = WCS(hsc_data['g']['hdr']).wcs.cd[1, 1] * 3600

    images = []
    # psfs = []  # TODO: include psf when available
    dp0_stats = []
    for band in 'grizy':
        images.append(hsc_data[band]['image'])
        # psfs.append(hsc_data[band]['psf'])  # TODO: include psf when available
        dp0_stats.append(dp0_sampler[band].sample()[0][0])

    # change zero points
    dp0_rms = [dp0_stats[b][1] for b in range(5)]
    dp0_zero_points = [dp0_stats[b][3] for b in range(5)]
    images = zero_point_change(images, 27, dp0_zero_points, dp0_rms, rms_frac_thresh=zp_rms_frac_thresh)

    degraded_images = []
    success = True
    for b, band in enumerate('grizy'):
        dp0_median, dp0_rms, dp0_seeing, dp0_zp = dp0_stats[b]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                deg_img_b = hsc_to_lsst(images[b],
                                        band,
                                        num_visits=num_visits,
                                        exp_time=30,
                                        lsst_zero_point=dp0_zp,
                                        hsc_fwhm=HSC_FWHM[band],
                                        # hsc_psf=psfs[b],  # TODO: include psf when available
                                        hsc_pix_scale=pix_scale,
                                        lsst_fwhm=dp0_seeing,
                                        background_noise=dp0_rms,
                                        background_median=dp0_median,
                                        psf_transform=True,
                                        add_poisson_noise=True,
                                        add_background_noise=True,
                                        use_nise_diff=True,
                                        to_adu=False,
                                        out_size=64
                                        )
            except Exception as e:
                print(f"Error in band {band}: {type(e)} {e}")
                success = False
                break
        degraded_images.append(deg_img_b)
    return success, degraded_images, dp0_zero_points


def write_degraded_image(
        filename,
        images,
        zero_points,
        exp_time=30,
        metadata=None
):
    images_stack = np.stack(images, axis=0)
    hdr = fits.Header()
    hdr['NAXIS'] = 3
    hdr['NAXIS1'] = images_stack.shape[1]
    hdr['NAXIS2'] = images_stack.shape[2]
    hdr['NAXIS3'] = len(images)
    hdr['BUNIT'] = 'cps'
    hdr['BAND_0'] = 'g'
    hdr['BAND_1'] = 'r'
    hdr['BAND_2'] = 'y'
    hdr['BAND_3'] = 'z'
    hdr['BAND_4'] = 'Y'
    hdr['ZP_G'] = zero_points[0]
    hdr['ZP_R'] = zero_points[1]
    hdr['ZP_I'] = zero_points[2]
    hdr['ZP_Z'] = zero_points[3]
    hdr['ZP_Y'] = zero_points[4]
    hdr['EXPTIME'] = exp_time
    if metadata is not None:
        for k, v in metadata.items():
            if isinstance(v, np.ma.core.MaskedConstant):
                v = None
            elif isinstance(v, str) and len(v) > 55:
                v = v[:55] + '...'
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                hdr[k.upper()] = v
    hdu = fits.PrimaryHDU(data=images_stack, header=hdr)
    hdu.writeto(filename, overwrite=True)


def query_degrade_write(
        ra,
        dec,
        size,
        field,
        username,
        password,
        dp0_sampler,
        out_filename,
        metadata=None,
        num_visits=1,
        zp_rms_frac_thresh=0.1
):
    degraded_images, success, zp = query_and_degrade(
        ra,
        dec,
        size,
        field,
        username,
        password,
        dp0_sampler,
        num_visits,
        zp_rms_frac_thresh
    )
    if success:
        write_degraded_image(
            out_filename,
            degraded_images,
            zp,
            exp_time=30 * num_visits,
            metadata=metadata
        )
    return success


