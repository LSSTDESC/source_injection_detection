import h5py
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.geom as geom
import os
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib.tools import *
import lsst.afw.image as afwImage

#============================
MAG0 = 28.17
LENS_FILENAME = "lens_data.h5"


#============================
def get_single_stamp(system_index, time_index, folder):

    system_tag = "system_%d"%system_index
    image_tag = "%s_%d"%(system_tag, time_index)

    with h5py.File(LENS_FILENAME, 'r') as h5f:
        
        system_metadata = pd.DataFrame(h5f[system_tag]["metadata"][:], 
                                       columns=h5f["metadata_columns"][:])

        # At that time point
        system_lc = h5f[system_tag]["light_curve"][:, time_index]
        system_image = h5f[system_tag]["images"][time_index]

    #----------------------------------
    total_flux = np.sum(system_image)
    total_mag = flux2mag(total_flux, MAG0)
    print("total_flux: ", total_flux)
    print("total_mag: ", total_mag)

    #----------------------------------
    lens_mag = system_metadata[b"deflector_light_magnitude"][0]
    lens_flux = mag2flux(lens_mag, MAG0)
    print("lens_flux: ", lens_flux)
    print("lens_mag: ", lens_mag)

    #----------------------------------
    source_number = len(system_lc)    
    source_flux_list = []

    for source_index in range(source_number):
        source_mag = system_lc[source_index]
        source_flux = mag2flux(source_mag, MAG0)
        source_flux_list.append(source_flux)
        print("source_flux: ", source_flux)
        print("source_mag: ", source_mag)
    
    #----------------------------------
    arc_flux = total_flux - lens_flux - np.sum(source_flux_list)
    arc_mag = flux2mag(arc_flux, MAG0)
    print("arc_flux: ", arc_flux)
    print("arc_mag: ", arc_mag)
    
    #----------------------------------
    plt.figure()
    #plt.imshow(np.log( system_image ), origin="lower")
    plt.imshow(np.arcsinh( system_image ), origin="lower")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("%s/%s.png"%(folder, image_tag) )
    plt.close()
    
    #----------------------------------
    fits.writeto("%s/%s.fits"%(folder, image_tag), system_image, overwrite=True)
    
    #----------------------------------
    
    return system_image


#============================
def get_coadd_stamp(system_index, folder):

    system_tag = "system_%d"%system_index
    image_tag = "%s_coadd"%system_tag

    with h5py.File(LENS_FILENAME, 'r') as h5f:
        system_image = h5f[system_tag]["images"][:]

    #----------------------------------
    mean_image = np.mean(system_image, axis=0)

    #----------------------------------
    total_flux = np.sum(mean_image)
    total_mag = flux2mag(total_flux, MAG0)
    print("total_flux: ", total_flux)
    print("total_mag: ", total_mag)

    #----------------------------------
    plt.figure()
    plt.imshow(np.arcsinh( mean_image ), origin="lower")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("%s/%s.png"%(folder, image_tag) )
    plt.close()
    
    #----------------------------------
    fits.writeto("%s/%s.fits"%(folder, image_tag), mean_image, overwrite=True)
    
    return mean_image


#============================
def get_diff_stamp(system_index, time_index, folder):

    soften_factor = 1.5

    system_tag = "system_%d"%system_index
    image_tag = "%s_diff"%system_tag

    #----------------------------------
    single_image_tag = "%s_%d"%(system_tag, time_index)
    single_image_path = "%s/%s.fits"%(folder, single_image_tag)
    
    if os.path.exists(single_image_path): 
        print("%s exists!"%single_image_path)
        single_image = fits.getdata(single_image_path)
        
    else:
        single_image = get_single_stamp(system_index, time_index)

    #----------------------------------
    coadd_image_tag = "%s_coadd"%system_tag
    coadd_image_path = "%s/%s.fits"%(folder, coadd_image_tag)
    
    if os.path.exists(coadd_image_path): 
        print("%s exists!"%coadd_image_path)
        coadd_image = fits.getdata(coadd_image_path)
        
    else:
        coadd_image = get_coadd_stamp(system_index)

    #----------------------------------
    diff_image = single_image - coadd_image
    threshold = np.max(np.abs(diff_image)) * soften_factor

    #----------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3), layout="tight")
    im = axs[0].imshow(diff_image, 
                       origin="lower", 
                       cmap="bwr", 
                       vmin=-threshold, 
                       vmax=threshold, 
                      )
    
    fig.colorbar(im, ax=axs[0])
    
    axs[0].set_title("Original")

    #----------------------------------
    sigma_pix = 2.
    pix2arcsec = 0.2
    # https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html
    seeing_pix = 2. * sigma_pix * np.sqrt( 2.*np.log(2.) )
    seeing_arcsec = seeing_pix * pix2arcsec # fwhm
    print("\nConvolving with seeing_arcsec: ", seeing_arcsec)

    diff_image_convolved = gaussian_filter(diff_image, sigma=sigma_pix)
    threshold = np.max(np.abs(diff_image_convolved)) * soften_factor
    
    im = axs[1].imshow(diff_image_convolved, 
                       origin="lower", 
                       cmap="bwr", 
                       vmin=-threshold, 
                       vmax=threshold, 
                      )
    
    fig.colorbar(im, ax=axs[1])

    axs[1].set_title("Conv %.2f\" Gaussian"%seeing_arcsec)
    
    #----------------------------------
    plt.savefig("%s/%s.png"%(folder, image_tag) )

    #----------------------------------

    return diff_image

    

#============================
def check_flux_diff(system_index):

    system_tag = "system_%d"%system_index

    with h5py.File(LENS_FILENAME, 'r') as h5f:
        system_lc = h5f[system_tag]["light_curve"][:, :]

    system_lc_flux = mag2flux(system_lc, MAG0)
    
    mean_flux = np.mean(system_lc_flux, axis=1)
    #print("mean_flux: ", mean_flux )

    diff_flux = system_lc_flux.T - mean_flux
    #print("diff_flux: \n", diff_flux )

    # There could be dipoles, so just consider sum
    sum_diff_flux = np.abs( np.sum(diff_flux, axis=1) )
    #print("\nsum_diff_flux: ", sum_diff_flux )

    argsort = np.argsort(sum_diff_flux)
    #print(argsort)
    max_index = argsort[-1]
    max_index2 = argsort[-2]
    print("max indices: ", max_index, max_index2)
    #print("max values: ", sum_diff_flux[max_index], sum_diff_flux[max_index2])
    
    return 0


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

def make_rotated_stamp(visit, stamp_filename):
    add_wcs(stamp_filename)
    wcs_stamp_filename = "%s_wcs.fits"%stamp_filename[:-5]
    stamp_img_orig = afwImage.ExposureF.readFits(wcs_stamp_filename)
    stamp_img_rotated = rotate_exposure(stamp_img_orig, -1*visit.visitInfo.getBoresightRotAngle().asDegrees())
    stamp_img_rotated.image.array[np.where(np.isnan(stamp_img_rotated.image.array))] = 0.0
    rot_stamp_filename = wcs_stamp_filename[:4]+'rotated_'+wcs_stamp_filename[4:]
    stamp_img_rotated.writeFits(rot_stamp_filename)