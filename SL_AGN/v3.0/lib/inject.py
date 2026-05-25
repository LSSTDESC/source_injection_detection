import lib.tools as tl
import lib.stamp as stp

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack

from lsst.source.injection import VisitInjectConfig, VisitInjectTask
from lsst.source.injection import CoaddInjectConfig, CoaddInjectTask



#======================================
def make_grid(ra_cen, dec_cen, width_arcmin, num_side):

    gap_arcsec = (width_arcmin * 60. - stp.STAMP_WIDTH_ARCSEC) / (num_side - 1) - stp.STAMP_WIDTH_ARCSEC
    print("gap_arcsec: ", gap_arcsec)
    
    if gap_arcsec < 6.0:
        print("ATT: Injection grid too tight!!")
        return None, None
    
    # We use plane sky approximation since it is a small region
    # But note angular length on RA needs a factor of cos(DEC)

    scale_factor = np.cos(np.deg2rad(dec_cen))
    ra_start = ra_cen - width_arcmin / 60. / 2 / scale_factor + 0.5 * stp.STAMP_WIDTH_ARCSEC / 3600.
    ra_end = ra_cen + width_arcmin / 60. / 2 / scale_factor - 0.5 * stp.STAMP_WIDTH_ARCSEC / 3600.

    dec_start = dec_cen - width_arcmin / 60. / 2 + 0.5 * stp.STAMP_WIDTH_ARCSEC / 3600.
    dec_end = dec_cen + width_arcmin / 60. / 2 - 0.5 * stp.STAMP_WIDTH_ARCSEC / 3600.

    ra_arr = np.linspace(ra_start, ra_end, num_side)
    dec_arr = np.linspace(dec_start, dec_end, num_side)

    RA_grid, DEC_grid = np.meshgrid(ra_arr, dec_arr)

    #offset = 1. / 3600. 
    #RA_grid += np.random.uniform(-offset, offset, (num_side, num_side)) / scale_factor
    #DEC_grid += np.random.uniform(-offset, offset, (num_side, num_side))

    RA_arr = RA_grid.flatten()
    DEC_arr = DEC_grid.flatten()

    #----------------------------------
    inj_radec = Table(
        {
            "ra": RA_arr,
            "dec": DEC_arr,
        }
    )
    
    #----------------------------------

    return inj_radec


def make_inj_catalog_stamp(wcs, stamp_mag_list, stamp_filename_list, inj_radec, tag):

    RA_arr, DEC_arr = inj_radec["ra"], inj_radec["dec"]
    x_arr, y_arr = wcs.skyToPixelArray(RA_arr, DEC_arr, degrees=True)

    n_inj = len(RA_arr)
    inj_catalog = Table(
        {
            "injection_id": range(n_inj),
            "ra": RA_arr,
            "dec": DEC_arr,
            "source_type": ["Stamp"] * n_inj,
            "mag": stamp_mag_list,
            "stamp": stamp_filename_list,
            "x": x_arr,
            "y": y_arr,
        }
    )

    #----------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 4),layout="constrained")
    
    axs[0].scatter(x_arr, y_arr, s=2.)
    
    # Draw a line pointing north
    ra_cen = np.median(RA_arr)
    dec_cen = np.median(DEC_arr)
    arrow_len_deg = 1. / 60
    dec_north = dec_cen + arrow_len_deg
    x_line_arr, y_line_arr = wcs.skyToPixelArray(np.array([ra_cen, ra_cen]),
                                                 np.array([dec_cen, dec_north]),
                                                 degrees=True)
    axs[0].annotate("", 
                 xytext=(x_line_arr[0], y_line_arr[0]), 
                 xy=(x_line_arr[1], y_line_arr[1]), 
                 arrowprops=dict(arrowstyle="->", color='red', linewidth=2))
    
    axs[0].set_xlabel("X [pix]")
    axs[0].set_ylabel("Y [pix]")

    axs[1].scatter(RA_arr, DEC_arr, s=1.)
    axs[1].set_xlabel("RA [deg]")
    axs[1].set_ylabel("DEC [deg]")
    axs[1].invert_xaxis()

    plt.savefig("%s/inj_%s.png"%(tl.FIG_FOLDER, tag) )

    #----------------------------------

    return inj_catalog 


def make_inj_catalog_visit(visit_image, stamp_mag_list, stamp_filename_list, inj_radec):

    wcs = visit_image.getWcs()

    md_dict = visit_image.metadata.toDict()
    visit = md_dict["LSST BUTLER DATAID VISIT"]
    detector = md_dict["LSST BUTLER DATAID DETECTOR"]

    visit_tag = "%d_%d"%(visit, detector)
    tag = "visit_%s"%visit_tag
    inj_catalog_visit = make_inj_catalog_stamp(wcs, stamp_mag_list, stamp_filename_list, inj_radec, tag)
    
    return inj_catalog_visit
    

def make_inj_catalog_template(template_image, stamp_mag_list, stamp_filename_list, inj_radec):

    wcs = template_image.getWcs()
    
    md_dict = template_image.getMetadata().toDict()
    tract = md_dict["LSST BUTLER DATAID TRACT"]
    patch = md_dict["LSST BUTLER DATAID PATCH"]
    band = md_dict["LSST BUTLER DATAID BAND"]
    
    template_tag = "%d_%d_%s"%(tract, patch, band)
    tag = "template_%s"%template_tag
    inj_catalog_template = make_inj_catalog_stamp(wcs, stamp_mag_list, stamp_filename_list, inj_radec, tag)
    
    return inj_catalog_template
    


#======================================
def image_inject(image, inj_catalog, image_type=None):

    if image_type=="visit":
        inject_config = VisitInjectConfig()
        inject_task = VisitInjectTask(config=inject_config)
    elif image_type=="template":
        inject_config = CoaddInjectConfig()
        inject_task = CoaddInjectTask(config=inject_config)
    else:
        print("ATT: Wrong image_type!")
        return None

    psf = image.getPsf()
    photo_calib = image.getPhotoCalib()
    wcs = image.getWcs()
    
    injected_output = inject_task.run(
        injection_catalogs=inj_catalog,
        input_exposure=image.clone(),
        psf=psf,
        photo_calib=photo_calib,
        wcs=wcs,
    )

    return injected_output.output_exposure, injected_output.output_catalog


#--------------------------------------
def visit_inject(visit_image, inj_catalog):

    return image_inject(visit_image, inj_catalog, image_type="visit")


def template_inject(template_image, inj_catalog):

    return image_inject(template_image, inj_catalog, image_type="template")
    


#======================================

def compute_coord_after_rot(x_c, y_c, delta_x, delta_y, rotation_angle):
    """Rotate (delta_x, delta_y) clockwise by rotation_angle (deg), then add to (x_c, y_c)."""
    theta = np.deg2rad(rotation_angle)
    delta_x_new = delta_x * np.cos(theta) + delta_y * np.sin(theta)
    delta_y_new = -delta_x * np.sin(theta) + delta_y * np.cos(theta)
    x = x_c + delta_x_new
    y = y_c + delta_y_new
    return x, y


def make_inj_catalog_point(wcs, inj_radec, point_mag_list, point_delta_x_list, point_delta_y_list, rotation_angle, tag, time_index_list=None):
    """
        point_mag_list, point_delta_x_list, point_delta_y_list: [[lens_i_image_0, lens_i_image_1, ...], [lens_i+1_image_0, lens_i+1_image_1, ...], ...]
        time_index_list: [time_index_for_lens_i, time_index_for_lens_i+1, ...] (one per system, -1 for template)
    """

    RA_arr = inj_radec["ra"]
    DEC_arr = inj_radec["dec"]
    x_c_arr, y_c_arr = wcs.skyToPixelArray(RA_arr, DEC_arr, degrees=True)

    all_ra = []
    all_dec = []
    all_mag = []
    all_x = []
    all_y = []
    all_time_index = []

    for i in range(len(inj_radec)):
        ti = time_index_list[i] if time_index_list is not None else -1
        for j in range(len(point_mag_list[i])):
            x, y = compute_coord_after_rot(x_c_arr[i], y_c_arr[i],
                                           point_delta_x_list[i][j],
                                           point_delta_y_list[i][j],
                                           rotation_angle)
            ra, dec = wcs.pixelToSkyArray(np.array([x]), np.array([y]), degrees=True)
            all_x.append(x)
            all_y.append(y)
            all_ra.append(ra[0])
            all_dec.append(dec[0])
            all_mag.append(point_mag_list[i][j])
            all_time_index.append(ti)

    n_point = len(all_ra)
    inj_catalog_point = Table(
        {
            "injection_id": range(n_point),
            "ra": all_ra,
            "dec": all_dec,
            "source_type": ["Star"] * n_point,
            "mag": all_mag,
            "x": all_x,
            "y": all_y,
            "time_index": all_time_index,
        }
    )

    return inj_catalog_point


def make_inj_catalog(wcs, stamp_mag_list, stamp_filename_list, inj_radec,
                     point_mag_list, point_delta_x_list, point_delta_y_list,
                     rotation_angle, tag, time_index_list=None):

    stamp_catalog = make_inj_catalog_stamp(wcs, stamp_mag_list, stamp_filename_list, inj_radec, tag)
    point_catalog = make_inj_catalog_point(wcs, inj_radec, point_mag_list,
                                           point_delta_x_list, point_delta_y_list,
                                           rotation_angle, tag,
                                           time_index_list=time_index_list)

    # Offset point source injection_ids to avoid overlap with stamp ids
    n_stamp = len(stamp_catalog)
    point_catalog["injection_id"] = point_catalog["injection_id"] + n_stamp

    inj_catalog = vstack([stamp_catalog, point_catalog])

    inj_catalog.write("%s/inj_catalog_%s.fits"%(tl.CATALOG_FOLDER, tag),
                      overwrite=True)

    return inj_catalog
