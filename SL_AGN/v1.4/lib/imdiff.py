from lsst.ip.diffim import subtractImages
from lsst.ip.diffim import detectAndMeasure
from lsst.ap.association import TransformDiaSourceCatalogTask
from lsst.pipe.base import Struct

import lsst.afw.geom
import lsst.meas.algorithms
import lsst.afw.math

import lib.tools as tl



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



def get_kernel(kernel_filename):
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

    initInputs = {}
    schema = diaSourceCat.getSchema()
    initInputs['diaSourceSchema'] = Struct(schema=schema)
    task = TransformDiaSourceCatalogTask(initInputs)
    result = task.run(diaSourceCat, diffIm, band, reliability)
    
    return result.diaSourceTable