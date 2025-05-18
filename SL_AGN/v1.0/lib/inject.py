from lib.tools import *

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table 

from lsst.source.injection import VisitInjectConfig, VisitInjectTask
from lsst.source.injection import CoaddInjectConfig, CoaddInjectTask



#======================================
SIZE_PIX = 4000


#======================================
def make_grid_coord(wcs, mag, stamp_filename, edge=200, step=500): 

    x = np.arange(edge, SIZE_PIX, step)
    y = np.arange(edge, SIZE_PIX, step)
    X, Y = np.meshgrid(x,y)

    n_inj = np.shape(X)[0] * np.shape(X)[1]

    x_list = X.flatten()
    y_list = Y.flatten()
  
    #----------------------------------
    ra_list = []
    dec_list = []

    for i in range(len(x_list)):
        x = x_list[i]
        y = y_list[i]

        tmp = wcs.pixelToSky(x, y)
        ra = tmp[0].asDegrees()
        dec = tmp[1].asDegrees()

        ra_list.append(ra)
        dec_list.append(dec)

    inj_catalog = Table(
        {
            "injection_id": range(n_inj),
            "ra": ra_list,
            "dec": dec_list,
            "source_type": ["Stamp"] * n_inj,
            "mag": [mag] * n_inj,
            "stamp": [stamp_filename] * n_inj,
            "x": x_list,
            "y": y_list,
        }
    )

    # May save the catalog to a file
    inj_catalog.write("%s/inj_catalog.fits"%CATALOG_FOLDER, 
                      overwrite=True)
    
    #----------------------------------
    fig, axs = plt.subplots(1,2,figsize=(8.5,4),layout="tight")
    axs[0].scatter(x_list, y_list)
    axs[0].set_xlabel("X [pix]")
    axs[0].set_ylabel("Y [pix]")
    axs[1].scatter(ra_list, dec_list)
    axs[1].set_xlabel("RA [deg]")
    axs[1].set_ylabel("DEC [deg]")
    axs[1].invert_xaxis()
    plt.savefig("%s/inj_coord.png"%FIG_FOLDER)

    #----------------------------------

    return inj_catalog 


#======================================
def calexp_inject_stamp(calexp, inj_catalog):

    psf = calexp.getPsf()
    photo_calib = calexp.getPhotoCalib()
    wcs = calexp.getWcs()

    inject_config = VisitInjectConfig()
    inject_task = VisitInjectTask(config=inject_config)

    injected_output = inject_task.run(
        injection_catalogs=inj_catalog,
        input_exposure=calexp.clone(),
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
