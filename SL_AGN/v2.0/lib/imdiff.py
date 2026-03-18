from lsst.ip.diffim import subtractImages
from lsst.ip.diffim import detectAndMeasure
from lsst.ap.association import TransformDiaSourceCatalogTask, TransformDiaSourceCatalogConfig
from lsst.pipe.base import Struct

import lsst.afw.geom
import lsst.meas.algorithms
import lsst.afw.math

import lib.tools as tl
import lib.visual as vis

from astropy.table import Table



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
    kernel = lsst.afw.math.LinearCombinationKernel.readFits(kernel_filename)
    return kernel

    

#=============================
def warp_image(visit_image, template_image):
    """
        Make a warped template.
    """

    config = lsst.afw.math.Warper.ConfigClass()
    config.warpingKernelName = "lanczos5"
    warper = lsst.afw.math.Warper.fromConfig(config)

    box = visit_image.getBBox()
    xyTransform = lsst.afw.geom.makeWcsPairTransform(template_image.wcs, visit_image.wcs)
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
def imdiff_detect(tract, patch, band, visit, detector, kernel=None):

    if kernel is None:
        prefix=''
    else:
        prefix='injected_'
    
    template_image_filename = f"{tl.FIG_FOLDER}/{prefix}template_{tract}_{patch}_{band}.fits"
    template_image = lsst.afw.image._exposure.ExposureF.readFits(template_image_filename)

    visit_image_filename = f"{tl.FIG_FOLDER}/{prefix}visit_{visit}_{detector}.fits"
    visit_image = lsst.afw.image._exposure.ExposureF.readFits(visit_image_filename)

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
    
    
    return kernel_new