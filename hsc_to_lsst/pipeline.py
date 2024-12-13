from hsc_to_lsst.hsc_query import query_hsc
from hsc_to_lsst.data_degradation.zero_point import zero_point_change
from hsc_to_lsst.data_degradation.hsc_degradation import hsc_to_lsst
from astropy.wcs import WCS
from astropy.io import fits
import warnings
# import tarfile


HSC_FWHM = {
    'g': 0.71,
    'r': 0.74,
    'i': 0.53,
    'z': 0.56,
    'y': 0.53
}


def query_and_degrade(
        ra,
        dec,
        dp0_sampler,
        semaphore,
        username=None,
        password=None,
        zp_rms_frac_thresh=0.3,
        hsc_size_arcsec=20,
        lsst_size_pix=61,
        field="pdr3_wide",
        verbose=False
):
    try:
        hsc_data = query_hsc(ra, dec, semaphore, username, password, hsc_size_arcsec, field)
    # this needs to keep going regardless of the error (common: OSError and tarfile.ReadError)
    except Exception as e:
        if verbose:
            print(f"Error querying HSC data: {type(e)} {e}")
        return False, None, None
    for band in 'grizy':
        if not hsc_data[band]:
            if verbose:
                print(f"Error querying HSC data: missing band {band}")
            return False, None, None

    pix_scale = WCS(hsc_data['g']['hdr']).wcs.cd[1, 1] * 3600

    images = []
    # psfs = []  # TODO: include psf when available
    dp0_stats = []
    for band in 'grizy':
        images.append(hsc_data[band]['image'])
        # psfs.append(hsc_data[band]['psf'])  # TODO: include psf when available
        dp0_stats.append(dp0_sampler[band].sample())

    # change zero points
    dp0_rms = [dp0_stats[b]['rms'] for b in range(5)]
    dp0_zero_points = [dp0_stats[b]['zero_point'] for b in range(5)]
    images, mag_change = zero_point_change(images,
                                           27,
                                           dp0_zero_points,
                                           dp0_rms,
                                           rms_frac_thresh=zp_rms_frac_thresh)

    degraded_images = []
    success = True
    for b, band in enumerate('grizy'):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                deg_img_b = hsc_to_lsst(images[b],
                                        band,
                                        exp_time=dp0_stats[b]['exp_time'],
                                        lsst_zero_point=dp0_stats[b]['zero_point'],
                                        hsc_fwhm=HSC_FWHM[band],
                                        # hsc_psf=psfs[b],  # TODO: include psf when available
                                        hsc_pix_scale=pix_scale,
                                        lsst_fwhm=dp0_stats[b]['seeing'],
                                        background_noise=dp0_stats[b]['rms'],
                                        background_median=dp0_stats[b]['median'],
                                        psf_transform=True,
                                        add_poisson_noise=True,
                                        add_background_noise=True,
                                        use_nise_diff=True,
                                        to_adu=False,
                                        out_size=lsst_size_pix
                                        )
            # this needs to keep going regardless of the error (common: ValueError)
            except Exception as e:
                if verbose:
                    print(f"Error processing band {band}: {type(e)} {e}")
                success = False
                break
        degraded_images.append(deg_img_b)
    return success, degraded_images, mag_change


def write_degraded_image(
        filename,
        images
):
    for b, band in enumerate('grizy'):
        hdu = fits.PrimaryHDU(data=images[b])
        hdu.writeto(f"{filename}_{band}.fits", overwrite=True)


def query_degrade_write(
        out_filename,
        ra,
        dec,
        dp0_sampler,
        semaphore,
        username=None,
        password=None,
        zp_rms_frac_thresh=0.3,
        hsc_size_arcsec=20,
        lsst_size_pix=61,
        field="pdr3_wide",
        verbose=False
):
    success, degraded_images, mag_change = query_and_degrade(
        ra,
        dec,
        dp0_sampler,
        semaphore,
        username,
        password,
        zp_rms_frac_thresh,
        hsc_size_arcsec,
        lsst_size_pix,
        field,
        verbose
    )
    if success:
        write_degraded_image(
            out_filename,
            degraded_images,
        )
    return success, mag_change
