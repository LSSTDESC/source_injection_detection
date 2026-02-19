import lib.tools as tl

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table 

from lsst.source.injection import VisitInjectConfig, VisitInjectTask
from lsst.source.injection import CoaddInjectConfig, CoaddInjectTask



#======================================
#SIZE_PIX = 4000
STAMP_WIDTH_PIX = 33
STAMP_WIDTH_ARCSEC = STAMP_WIDTH_PIX / tl.ARCSEC2PIX


#======================================
def make_grid(ra_cen, dec_cen, width_arcmin, num_side):

    gap_arcsec = (width_arcmin * 60. - STAMP_WIDTH_ARCSEC) / (num_side - 1) - STAMP_WIDTH_ARCSEC
    
    #test_side_length = (stamp_width_pix + gap_pix) * (num_side - 1) + stamp_width_pix
    
    #gap_pix = gap_arcsec * tl.ARCSEC2PIX
    
    #if gap_pix < (0.5 * stamp_width_pix):
    
    if gap_arcsec < 6.0:
        print("ATT: Injection grid too tight!!")
        return None, None
    
    # We use plane sky approximation since it is a small region
    # But note angular length on RA needs a factor of cos(DEC)

    ra_start = ra_cen - width_arcmin / 60. / 2 / np.cos(np.deg2rad(dec_cen)) + 0.5 * STAMP_WIDTH_ARCSEC / 3600.
    ra_end = ra_cen + width_arcmin / 60. / 2 / np.cos(np.deg2rad(dec_cen)) - 0.5 * STAMP_WIDTH_ARCSEC / 3600.

    dec_start = dec_cen - width_arcmin / 60. / 2 + 0.5 * STAMP_WIDTH_ARCSEC / 3600.
    dec_end = dec_cen + width_arcmin / 60. / 2 - 0.5 * STAMP_WIDTH_ARCSEC / 3600.

    ra_arr = np.linspace(ra_start, ra_end, num_side)
    dec_arr = np.linspace(dec_start, dec_end, num_side)

    RA_grid, DEC_grid = np.meshgrid(ra_arr, dec_arr)

    RA_arr = RA_grid.flatten()
    DEC_arr = DEC_grid.flatten()

    #----------------------------------
    inj_radec = Table(
        {
            "ra": RA_arr,
            "dec": DEC_arr,
        }
    )
    inj_radec.write("%s/inj_radec.fits"%tl.FIG_FOLDER, overwrite=True)
    
    #----------------------------------
    fig, axs = plt.subplots(1, 1, figsize=(4, 4), layout="constrained")
    axs.scatter(RA_arr, DEC_arr)
    axs.set_xlabel("RA [deg]")
    axs.set_ylabel("DEC [deg]")
    axs.invert_xaxis()
    plt.savefig("%s/inj_radec.png"%tl.FIG_FOLDER)

    #----------------------------------

    return RA_arr, DEC_arr


def make_inj_catalog(wcs, stamp_mag, stamp_filename, RA_arr, DEC_arr, tag):

    x_arr, y_arr = wcs.skyToPixelArray(RA_arr, DEC_arr, degrees=True)

    n_inj = len(RA_arr)
    inj_catalog = Table(
        {
            "injection_id": range(n_inj),
            "ra": RA_arr,
            "dec": DEC_arr,
            "source_type": ["Stamp"] * n_inj,
            "mag": [stamp_mag] * n_inj,
            "stamp": [stamp_filename] * n_inj,
            #"x": x_arr,
            #"y": y_arr,
        }
    )

    #----------------------------------
    inj_catalog.write("%s/inj_catalog_%s.fits"%(tl.FIG_FOLDER, tag), 
                      overwrite=True)
    
    #----------------------------------
    fig, axs = plt.subplots(1, 1, figsize=(4, 4),layout="constrained")
    
    axs.scatter(x_arr, y_arr)
    
    # Draw a line pointing north
    ra_cen = np.median(RA_arr)
    dec_cen = np.median(DEC_arr)
    arrow_len_deg = 0.2 / 60
    dec_north = dec_cen + arrow_len_deg
    x_line_arr, y_line_arr = wcs.skyToPixelArray(np.array([ra_cen, ra_cen]),
                                                 np.array([dec_cen, dec_north]),
                                                 degrees=True)
    axs.annotate("", 
                 xytext=(x_line_arr[0], y_line_arr[0]), 
                 xy=(x_line_arr[1], y_line_arr[1]), 
                 arrowprops=dict(arrowstyle="->", color='red', linewidth=2))
    
    axs.set_xlabel("X [pix]")
    axs.set_ylabel("Y [pix]")

    plt.savefig("%s/inj_xy_%s.png"%(tl.FIG_FOLDER, tag) )

    #----------------------------------

    return inj_catalog, x_arr, y_arr 


def make_inj_catalog_visit(visit_image, stamp_mag, stamp_filename, RA_arr, DEC_arr):

    wcs = visit_image.getWcs()

    md_dict = visit_image.metadata.toDict()
    visit = md_dict["LSST BUTLER DATAID VISIT"]

    visit_tag = "%d"%visit
    tag = "visit_%s"%visit_tag
    inj_catalog_visit = make_inj_catalog(wcs, stamp_mag, stamp_filename, RA_arr, DEC_arr, tag)
    
    return inj_catalog_visit
    

def make_inj_catalog_template(template_image, stamp_mag, stamp_filename, RA_arr, DEC_arr):

    wcs = template_image.getWcs()
    
    md_dict = template_image.getMetadata().toDict()
    tract = md_dict["LSST BUTLER DATAID TRACT"]
    patch = md_dict["LSST BUTLER DATAID PATCH"]
    
    template_tag = "%d_%d"%(tract, patch)
    tag = "template_%s"%template_tag
    inj_catalog_template = make_inj_catalog(wcs, stamp_mag, stamp_filename, RA_arr, DEC_arr, tag)
    
    return inj_catalog_template
    


#======================================
def image_inject_stamp(image, inj_catalog, image_type=None):

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
    #print(injected_output.output_catalog)
    #injected_exposure = injected_output.output_exposure

    return injected_output.output_exposure, injected_output.output_catalog


def visit_inject_stamp(visit_image, inj_catalog):

    return image_inject_stamp(visit_image, inj_catalog, image_type="visit")


#--------------------------------------
def template_inject_stamp(template_image, inj_catalog):

    return image_inject_stamp(template_image, inj_catalog, image_type="template")
    
#-------------------
