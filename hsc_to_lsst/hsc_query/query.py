from hsc_to_lsst.hsc_query import downloadCutout, downloadPsf
import io
from astropy.io import fits


def query_hsc(username, password, ra, dec, size=10, field="pdr3_wide"):
    rect = downloadCutout.Rect.create(
        ra=str(ra),
        dec=str(dec),
        sw=f"{size/2}arcmin",
        sh=f"{size/2}arcmin",
        filter="all",
        rerun=field
    )
    psf_req = downloadPsf.PsfRequest.create(
        ra=str(ra),
        dec=str(dec),
        filter="all",
        rerun=field,
    )
    image_list = downloadCutout.download(rect, user=username, password=password)
    # psf_list = downloadPsf.download(psf_req, user=username, password=password)

    output_data = {}
    for band in 'grizy':
        band_data = {}
        for metadata, data in image_list:
            if metadata['filter'] == f"HSC-{band.upper()}":
                band_data['image'], band_data['hdr'] = fits.getdata(io.BytesIO(data), header=True)
                band_data['metadata'] = metadata
                break
        # for psf_metadata, psf_data in psf_list:
        #     if psf_metadata['filter'] == f"HSC-{filter_name.upper()}":
        #         psf = fits.getdata(io.BytesIO(psf_data))
        #         filter_data['psf'] = psf
        #         break
        output_data[band] = band_data
    return output_data
