from lsst.ip.diffim import subtractImages
from lsst.ip.diffim import detectAndMeasure
from lsst.ap.association import TransformDiaSourceCatalogTask, TransformDiaSourceCatalogConfig
from lsst.pipe.base import Struct

import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.meas.algorithms

import lib.tools as tl
import lib.visual as vis

from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS as astropyWCS
import numpy as np



#=============================
def make_subtract_task(useExistingKernel=False):
    """
        Set up task for image subtraction between original images (no existing kernel).

        https://github.com/lsst/ip_diffim/blob/main/tests/test_subtractTask.py#L1339
    """

    subtractTask = subtractImages.SimplifiedSubtractTask

    #fluxField = "truth_instFlux"
    #errField = "truth_instFluxErr"

    fluxField = "base_PsfFlux_instFlux"
    errField = "base_PsfFlux_instFluxErr"

    config = subtractTask.ConfigClass()
    config.doSubtractBackground = False
    config.restrictKernelEdgeSources = False
    config.sourceSelector.signalToNoise.fluxField = fluxField
    config.sourceSelector.signalToNoise.errField = errField
    config.sourceSelector.doUnresolved = True
    config.sourceSelector.doIsolated = True
    config.sourceSelector.doRequirePrimary = True
    config.sourceSelector.doFlags = True
    config.sourceSelector.doSignalToNoise = True
    config.sourceSelector.flags.bad = ["base_PsfFlux_flag", ]

    config.update(useExistingKernel=useExistingKernel)
    task = subtractTask(config=config)

    return task



def run_subtract_task(warped_template, visit_image, kernel=None):

    if kernel is None:
        useExistingKernel = False
        task = make_subtract_task(useExistingKernel=useExistingKernel)
        result = task.run(warped_template, visit_image)
    else:
        useExistingKernel = True
        task = make_subtract_task(useExistingKernel=useExistingKernel)
        result = task.run(warped_template, visit_image, inputPsfMatchingKernel=kernel)
        
    difference_image = result.difference
    kernel = result.psfMatchingKernel

    return difference_image, kernel



def read_kernel(kernel_filename):
    kernel = afwMath.LinearCombinationKernel.readFits(kernel_filename)
    return kernel

    

#=============================
def warp_image(visit_image, template_image):
    """
        Make a warped template.
    """

    config = afwMath.Warper.ConfigClass()
    config.warpingKernelName = "lanczos5"
    warper = afwMath.Warper.fromConfig(config)

    box = visit_image.getBBox()
    xyTransform = afwGeom.makeWcsPairTransform(template_image.wcs, visit_image.wcs)
    warpedPsf = lsst.meas.algorithms.WarpedPsf(template_image.psf, xyTransform)
    warped = warper.warpExposure(visit_image.wcs, template_image, destBBox=box)
    warped.setPsf(warpedPsf)

    return warped



#=============================
def detect_dia_sources(visit_image, warped_template, difference_image): 

    task = detectAndMeasure.DetectAndMeasureTask()
    
    result = task.run(visit_image, warped_template, difference_image)

    dia_sources = result.diaSources #.asAstropy()

    return dia_sources


def sourceCat2Tab(diaSourceCat, diffIm, band, reliability=None):

    """
        Following the test. 
        https://github.com/lsst/ap_association/blob/main/tests/test_transformDiaSourceCatalog.py#L95
    """

    initInputs = {}
    schema = diaSourceCat.getSchema()
    
    initInputs['diaSourceSchema'] = Struct(schema=schema)
    config = TransformDiaSourceCatalogConfig()
    #config.doIncludeReliability = True
    task = TransformDiaSourceCatalogTask(initInputs, config=config)
    
    result = task.run(diaSourceCat, diffIm, band) #, reliability)
    
    return result.diaSourceTable


#============================
def _get_cutout_bounds(x, y, half, bbox):
    """Get clamped cutout bounds within image bounding box."""
    xmin = max(int(x - half), bbox.beginX)
    xmax = min(int(x + half), bbox.endX)
    ymin = max(int(y - half), bbox.beginY)
    ymax = min(int(y + half), bbox.endY)
    return xmin, xmax, ymin, ymax


def save_cutouts_individual(image, inj_radec, tag):
    """Save one FITS file per injected source (image + variance + mask)."""

    wcs = image.getWcs()
    RA_arr = np.array(inj_radec['ra'])
    DEC_arr = np.array(inj_radec['dec'])
    x_arr, y_arr = wcs.skyToPixelArray(RA_arr, DEC_arr, degrees=True)

    half = tl.CUTOUT_SIZE // 2
    bbox = image.getBBox()

    for i in range(len(inj_radec)):
        xmin, xmax, ymin, ymax = _get_cutout_bounds(x_arr[i], y_arr[i], half, bbox)
        cutout = image[xmin:xmax, ymin:ymax]
        cutout.writeFits(f"{tl.CUTOUT_FOLDER}/{tag}_id{i}.fits")


def save_cutouts_stacked(image, inj_radec, tag):
    """Save 3 multi-extension FITS files (image, variance, mask), one extension per system."""

    wcs = image.getWcs()
    RA_arr = np.array(inj_radec['ra'])
    DEC_arr = np.array(inj_radec['dec'])
    x_arr, y_arr = wcs.skyToPixelArray(RA_arr, DEC_arr, degrees=True)

    half = tl.CUTOUT_SIZE // 2
    bbox = image.getBBox()

    # Get parent WCS as astropy header, and parent xy0 for offset calculation
    wcs_header = astropyWCS(image.getWcs().getFitsMetadata()).to_header()
    image_xy0 = image.getXY0()

    hdus_image = [fits.PrimaryHDU()]
    hdus_variance = [fits.PrimaryHDU()]
    hdus_mask = [fits.PrimaryHDU()]

    for i in range(len(inj_radec)):
        xmin, xmax, ymin, ymax = _get_cutout_bounds(x_arr[i], y_arr[i], half, bbox)
        cutout = image[xmin:xmax, ymin:ymax]

        xy0 = cutout.getXY0()

        # Adjust CRPIX from parent image frame to cutout LOCAL frame
        cutout_wcs_header = wcs_header.copy()
        cutout_wcs_header['CRPIX1'] -= (xy0.getX() - image_xy0.getX())
        cutout_wcs_header['CRPIX2'] -= (xy0.getY() - image_xy0.getY())

        header_kwargs = {
            'INJ_ID': (i, 'injection id'),
            'RA': (RA_arr[i], 'RA in degrees'),
            'DEC': (DEC_arr[i], 'DEC in degrees'),
            'X0': (xy0.getX(), 'PARENT x origin (xy0)'),
            'Y0': (xy0.getY(), 'PARENT y origin (xy0)'),
        }

        hdu_img = fits.ImageHDU(cutout.image.array, name=f"ID_{i:04d}")
        hdu_var = fits.ImageHDU(cutout.variance.array, name=f"ID_{i:04d}")
        hdu_msk = fits.ImageHDU(cutout.mask.array, name=f"ID_{i:04d}")

        for hdu in [hdu_img, hdu_var, hdu_msk]:
            hdu.header.update(cutout_wcs_header)
            for key, val in header_kwargs.items():
                hdu.header[key] = val

        hdus_image.append(hdu_img)
        hdus_variance.append(hdu_var)
        hdus_mask.append(hdu_msk)

    fits.HDUList(hdus_image).writeto(f"{tl.CUTOUT_FOLDER}/{tag}_image.fits", overwrite=True)
    fits.HDUList(hdus_variance).writeto(f"{tl.CUTOUT_FOLDER}/{tag}_variance.fits", overwrite=True)
    fits.HDUList(hdus_mask).writeto(f"{tl.CUTOUT_FOLDER}/{tag}_mask.fits", overwrite=True)


#============================
def imdiff_detect(tract, patch, band, visit, detector, kernel=None,
                  if_save_cutout_individual=False, if_save_cutout_stacked=False):

    if kernel is None:
        prefix=''
    else:
        prefix='injected_'
    
    template_image_filename = f"{tl.FIG_FOLDER}/{prefix}template_{tract}_{patch}_{band}.fits"
    template_image = afwImage.ExposureF.readFits(template_image_filename)

    visit_image_filename = f"{tl.FIG_FOLDER}/{prefix}visit_{visit}_{detector}.fits"
    visit_image = afwImage.ExposureF.readFits(visit_image_filename)

    # Warp the template to match the visit image.
    warped_template = warp_image(visit_image, template_image)

    # Run image subtraction and save the kernel.
    difference_image, kernel_new = run_subtract_task(warped_template, visit_image, kernel)

    inj_radec = Table.read(f"{tl.CATALOG_FOLDER}/inj_catalog_visit_{visit}_{detector}.fits")
    
#    if prefix=='':
#        tag = "original"
#    elif prefix=='injected_':
#        tag = "injected"
#    else:
#        print("ATT: Wrong image_type!")
#        return 1

    tag = f"{visit}_{detector}_{tract}_{patch}_{band}"
    
    vis.plot_triple(visit_image, warped_template, difference_image, inj_radec, f"{prefix}{tag}")

    difference_image.writeFits(f"{tl.FIG_FOLDER}/{prefix}diff_image_{tag}.fits")

    dia_sources = detect_dia_sources(visit_image, warped_template, difference_image)

    dia_sources.writeFits(f"{tl.CATALOG_FOLDER}/{prefix}diaSources_{tag}.fits")

    dia_sources_tab_df = sourceCat2Tab(dia_sources, difference_image, difference_image.getFilter())
    dia_sources_tab_df.to_csv(f"{tl.CATALOG_FOLDER}/{prefix}diaSources_tab_{tag}.csv", index=False)

    # Save cutouts for injected images
    if kernel is not None:
        images_and_tags = [
            (visit_image, f"visit_{tag}"),
            (warped_template, f"template_{tag}"),
            (difference_image, f"diff_{tag}"),
        ]
        for img, cutout_tag in images_and_tags:
            if if_save_cutout_individual:
                save_cutouts_individual(img, inj_radec, cutout_tag)
            if if_save_cutout_stacked:
                save_cutouts_stacked(img, inj_radec, cutout_tag)

    return kernel_new