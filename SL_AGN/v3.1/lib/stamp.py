import lib.tools as tl
import h5py
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.geom as geom



#============================
# https://smtn-002.lsst.io/
MAG_ZPS_CPS = {
'u': 26.52,
'g': 28.51,
'r': 28.36,
'i': 28.17,
'z': 27.78,
'y': 26.82,
}


LENS_FILENAMES = [
    "3000sqdeg_lsst_1y_sample_merged_cleaned.h5",
]


STAMP_WIDTH_PIX = 33
STAMP_WIDTH_ARCSEC = STAMP_WIDTH_PIX / tl.ARCSEC2PIX


#============================
def get_n_observations(system_index, band, lens_filename, id_offset=0):
    """Return number of observation epochs, or None if system exceeds limit (>30)."""
    hdf5_index = system_index - id_offset
    with h5py.File(lens_filename, 'r') as h5f:
        try:
            n = len(h5f[f"lsst_lens_{hdf5_index}"]['observation_dates'][band])
        except KeyError:
            return None
    return n if n <= 30 else None


def get_lens_file_info():
    """Return list of (filename, id_offset, n_systems) for each lens file."""
    info = []
    id_offset = 0
    for filename in LENS_FILENAMES:
        with h5py.File(filename, 'r') as h5f:
            n_systems = len(h5f.keys())
        info.append((filename, id_offset, n_systems))
        id_offset += n_systems
    return info

#============================
def flux2mag(flux, band, exposure_time=30.):
    output = -2.5 * np.log10(flux/exposure_time) + MAG_ZPS_CPS[band]
    return output


def mag2flux(mag, band, exposure_time=30.):
    flux = 10**((mag - MAG_ZPS_CPS[band]) / (-2.5)) * exposure_time
    return flux



#======================================
def add_wcs(filename):

    data, hdr = fits.getdata(filename, ext=0, header=True)
    shape = np.shape(data)
    nrow = shape[0]
    ncol = shape[1]

    #----------------------------------
    w = WCS(naxis=2)
    w.wcs.crpix = [nrow/2., ncol/2.] # this is from shape of data
    w.wcs.crval = [0.0, 0.0]
    w.wcs.cdelt = np.array([-5.55555555555556E-05,5.55555555555556E-05])
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    wcs_header = w.to_header()

    #----------------------------------
    tag = filename.split('.')[-2]
    
    hdr["EXTNAME"] = "IMAGE"
    hdr_new = hdr + wcs_header

    #----------------------------------
    fits.writeto("%s_wcs.fits"%tag, data, hdr_new, overwrite=True)
    
    return 0

    
def rotate_exposure(exp, n_degrees):

    n_degrees = n_degrees % 360

    wcs = exp.getWcs()

    warper = afwMath.Warper('lanczos4')

    affine_rot_transform = geom.AffineTransform.makeRotation(n_degrees*geom.degrees)
    transform_p2top2 = afwGeom.makeTransform(affine_rot_transform)
    rotated_wcs = afwGeom.makeModifiedWcs(transform_p2top2, wcs, False)

    rotated_exp = warper.warpExposure(rotated_wcs, exp)

    return rotated_exp


def make_rotated_stamp(rotation_angle, wcs_stamp_filename, rot_wcs_stamp_filename):
    
    stamp_img_orig = afwImage.ExposureF.readFits(wcs_stamp_filename)
    
    stamp_img_rotated = rotate_exposure(stamp_img_orig, rotation_angle)
    
    stamp_img_rotated.image.array[np.where(np.isnan(stamp_img_rotated.image.array))] = 0.0
    
    stamp_img_rotated.writeFits(rot_wcs_stamp_filename)
    
    return 0
    


#======================================
def get_static_stamp(system_index, band, folder, lens_filename=None, id_offset=0):

    if lens_filename is None:
        lens_filename = LENS_FILENAMES[0]

    hdf5_index = system_index - id_offset

    system_tag = f"system_{system_index}"
    image_tag = f"{system_tag}_static_{band}"

    with h5py.File(lens_filename, 'r') as h5f:
        stamp_image = h5f[f"lsst_lens_{hdf5_index}"]['static_image'][band]['lens_plus_lensed_agn_host'][:]

    hdu = fits.PrimaryHDU(stamp_image)
    hdu.header["TOT_MAG"] = flux2mag(np.sum(stamp_image), band)
    hdu.writeto(f"{folder}/{image_tag}.fits", overwrite=True)

    return stamp_image


def get_point_mags(system_index, time_index, band, lens_filename=None, id_offset=0):

    if lens_filename is None:
        lens_filename = LENS_FILENAMES[0]

    hdf5_index = system_index - id_offset

    with h5py.File(lens_filename, 'r') as h5f:
        light_curves = h5f[f"lsst_lens_{hdf5_index}"]['light_curves']
        n_images = len(light_curves)
        mags = []
        for image_ind in range(n_images):
            mag = light_curves[f'image_{image_ind}'][band][time_index]
            mags.append(mag)

    return mags


def get_point_coords(system_index, band, lens_filename=None, id_offset=0):

    if lens_filename is None:
        lens_filename = LENS_FILENAMES[0]

    hdf5_index = system_index - id_offset

    with h5py.File(lens_filename, 'r') as h5f:
        metadata = h5f[f"lsst_lens_{hdf5_index}"]['metadata']
        light_curves = h5f[f"lsst_lens_{hdf5_index}"]['light_curves']
        n_images = len(light_curves)
        delta_x_list = []
        delta_y_list = []
        for ind in range(n_images):
            delta_x = metadata.attrs[f'point_source_light_{band}_ra_image_{ind}'] / tl.PIX2ARCSEC
            delta_y = metadata.attrs[f'point_source_light_{band}_dec_image_{ind}'] / tl.PIX2ARCSEC
            delta_x_list.append(delta_x)
            delta_y_list.append(delta_y)

    return delta_x_list, delta_y_list


#======================================
def get_all_point_delta_coords(system_indices, band, lens_filename=None, id_offset=0):

    all_delta_x = []
    all_delta_y = []
    for system_index in system_indices:
        delta_x, delta_y = get_point_coords(system_index, band,
                                            lens_filename=lens_filename,
                                            id_offset=id_offset)
        all_delta_x.append(delta_x)
        all_delta_y.append(delta_y)

    return all_delta_x, all_delta_y


def get_all_point_mags_visit(system_indices, time_indices, band, lens_filename=None, id_offset=0):

    all_mags = []
    for i, system_index in enumerate(system_indices):
        mags = get_point_mags(system_index, time_indices[i], band,
                              lens_filename=lens_filename,
                              id_offset=id_offset)
        all_mags.append(mags)

    return all_mags


def get_all_point_mags_template(system_indices, band, lens_filename=None, id_offset=0):

    if lens_filename is None:
        lens_filename = LENS_FILENAMES[0]

    all_mean_mags = []
    with h5py.File(lens_filename, 'r') as h5f:
        for system_index in system_indices:
            hdf5_index = system_index - id_offset
            light_curves = h5f[f"lsst_lens_{hdf5_index}"]['light_curves']
            n_images = len(light_curves)
            mean_mags = []
            for image_ind in range(n_images):
                mags_all_epochs = light_curves[f'image_{image_ind}'][band][:]
                fluxes = mag2flux(mags_all_epochs, band)
                mean_flux = np.mean(fluxes)
                mean_mag = flux2mag(mean_flux, band)
                mean_mags.append(mean_mag)
            all_mean_mags.append(mean_mags)

    return all_mean_mags
