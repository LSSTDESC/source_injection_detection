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

    ra_start = ra_cen - width_arcmin / 60. / 2 / np.cos(np.deg2rad(dec_center)) + 0.5 * STAMP_WIDTH_ARCSEC
    ra_end = ra_cen + width_arcmin / 60. / 2 / np.cos(np.deg2rad(dec_center)) - 0.5 * STAMP_WIDTH_ARCSEC

    dec_start = dec_cen - width_arcmin / 60. / 2 + 0.5 * STAMP_WIDTH_ARCSEC
    dec_end = dec_cen + width_arcmin / 60. / 2 - 0.5 * STAMP_WIDTH_ARCSEC

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
            "x": x_arr,
            "y": y_arr,
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
    dec_north = dec_cen + 0.5 / 60
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

    return inj_catalog 


def make_inj_catalog_visit(visit_wcs, stamp_mag, stamp_filename, RA_arr, DEC_arr, visit_tag):

    wcs = visit_wcs
    tag = "visit_%s"%visit_tag
    inj_catalog_visit = make_inj_catalog(wcs, stamp_mag, stamp_filename, RA_arr, DEC_arr, tag)
    
    return inj_catalog_visit
    

def make_inj_catalog_template(template_wcs, stamp_mag, stamp_filename, RA_arr, DEC_arr, template_tag):

    wcs = template_wcs
    tag = "template_%s"%visit_tag
    inj_catalog_visit = make_inj_catalog(wcs, stamp_mag, stamp_filename, RA_arr, DEC_arr, tag)
    
    return inj_catalog_template
    


#======================================

def visit_inject_stamp(visit, inj_catalog):

    psf = visit.getPsf()
    photo_calib = visit.getPhotoCalib()
    wcs = visit.getWcs()

    inject_config = VisitInjectConfig()
    inject_task = VisitInjectTask(config=inject_config)

    injected_output = inject_task.run(
        injection_catalogs=inj_catalog,
        input_exposure=visit.clone(),
        psf=psf,
        photo_calib=photo_calib,
        wcs=wcs,
    )
    injected_exposure = injected_output.output_exposure

    return injected_exposure


#--------------------------------------
def template_inject_stamp(template, inj_catalog):

    psf = template.getPsf()
    photo_calib = template.getPhotoCalib()
    wcs = template.getWcs()

    inject_config = CoaddInjectConfig()
    inject_task = CoaddInjectTask(config=inject_config)

    injected_output = inject_task.run(
        injection_catalogs=inj_catalog,
        input_exposure=template.clone(),
        psf=psf,
        photo_calib=photo_calib,
        wcs=wcs,
    )
    injected_exposure = injected_output.output_exposure

    return injected_exposure

#-------------------
