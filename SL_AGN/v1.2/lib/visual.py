from lib.tools import *

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.visualization import AsinhStretch, LinearStretch, LogStretch, ImageNormalize

import lsst.afw.display as afwDisplay



#============================
afwDisplay.setDefaultBackend("matplotlib")
#plt.style.use('tableau-colorblind10')


#============================
def clean_ax(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    return 0
    

def display_lsp(image, fig, ax):

    # Input: image from the LSST Science Pipelines

    #fig, ax = plt.subplots(figsize=(4,4))

    plt.sca(ax)

    display = afwDisplay.Display(frame=fig)
    display.scale("linear", "zscale")
    display.mtv(image.image)

    return display
    

def get_astropy_wcs(image):

    # Input: image from the LSST Science Pipelines

    w0 = image.getWcs()
    w = WCS(naxis=2)
    w.wcs.crpix = np.array(w0.getPixelOrigin()) # pix
    w.wcs.crval = [w0.getSkyOrigin()[0].asDegrees(),
                   w0.getSkyOrigin()[1].asDegrees()] # ra, dec
    w.wcs.cd = w0.getCdMatrix()
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    return w


def display_plt(image, fig, pos, stretch="log"):

    # Input: image from the LSST Science Pipelines

    #fig = plt.figure(figsize=(4,4))

    wcs = get_astropy_wcs(image)

    ax = fig.add_subplot(pos, projection=wcs)
    
    clean_ax(ax)

    fmin = np.min(image.image.array)
    fmax = np.max(image.image.array)
    fstd = np.std(image.image.array)
    fmed = np.median(image.image.array)
    print("min: ", fmin )
    print("max: ", fmax )
    print("std: ", fstd )
    print("med: ", fmed )
    vmax = fmed + fstd * 2.
    vmin = fmed - fstd * 1.

    if stretch=="asinh":
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(a=1e0))
    elif stretch=="linear":
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
    elif stretch=="log":
        norm = ImageNormalize(vmin=1e-3, vmax=vmax, stretch=LogStretch())
    else:
        print("ATT: stretch is wrong!")
        return 1
    
    im = ax.imshow(image.image.array, norm=norm) #, origin='lower')

    for c in ax.coords:
        c.set_major_formatter("d.dd")
        c.set_ticks(number=10)

    ax.grid(ls=':')
    #ax.set_xlabel('RA')
    #ax.set_ylabel('DEC')

    fig.colorbar(im, ax=ax)

    return 0


def plot_save_image(image, tag, method=None, plt_stretch="log"):
    
    if method is None: 
        
        fig = plt.figure(figsize=(10, 4), dpi=100)
        ax1 = fig.add_subplot(1, 2, 1)
        display_lsp(image, fig, ax1)
        
        display_plt(image, fig, 122)
        plt.tight_layout()
        plt.savefig("%s/%s_lp.png"%(FIG_FOLDER, tag))

    elif method=="lsp":
        
        fig, ax = plt.subplots(1, 1, dpi=100)
        display_lsp(image, fig, ax)
        plt.tight_layout()
        plt.savefig("%s/%s_lsp.png"%(FIG_FOLDER, tag))

    elif method=="plt":
        
        fig = plt.figure(dpi=100)
        display_plt(image, fig, 111, stretch=plt_stretch)
        plt.tight_layout()
        plt.savefig("%s/%s_plt.png"%(FIG_FOLDER, tag))

    else:
        print("ATT: wrong method!")
        return 1

    
    return 0


def plot_save_two_images(image1, image2, tag1, tag2, method="lsp"):
    
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=100)

    if method=="lsp":
        display_lsp(image1, fig, axs[0])
        axs[0].set_title(tag1)

        display_lsp(image2, fig, axs[1])
        axs[1].set_title(tag2)

    else:
        print("ATT: wrong method!")
        return 1

    #plt.tight_layout()
    plt.savefig("%s/%s-%s.png"%(FIG_FOLDER, tag1, tag2), bbox_inches="tight", pad_inches=0.1 )

    return 0


def plot_save_three_images(image1, image2, image3, tag1, tag2, tag3, method="lsp", xlim=None, ylim=None):
    
    fig, axs = plt.subplots(1, 3, figsize=(15, 4), layout="constrained")

    if method=="lsp":
        display_lsp(image1, fig, axs[0])
        axs[0].set_title(tag1)

        display_lsp(image2, fig, axs[1])
        axs[1].set_title(tag2)

        display_lsp(image3, fig, axs[2])
        axs[2].set_title(tag3)

        if xlim is not None and ylim is not None:
            for ax in axs:
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

        outname = "%s/%s-%s-%s_lsp.png"%(FIG_FOLDER, tag1, tag2, tag3)

    elif method=="plt" and xlim is not None and ylim is not None:
        xmin, xmax = xlim[0], xlim[1]
        ymin, ymax = ylim[0], ylim[1]
        im = axs[0].imshow(np.arcsinh(image1.image.array[ymin+10:ymax+10, xmin+10:xmax+10]), origin='lower')
        axs[0].set_title(tag1)
        fig.colorbar(im, ax=axs[0])#, orientation="horizontal")

        im = axs[1].imshow(np.arcsinh(image2.image.array[ymin:ymax, xmin:xmax]), origin='lower')
        axs[1].set_title(tag2)
        fig.colorbar(im, ax=axs[1])#, orientation="horizontal")

        im = axs[2].imshow(np.arcsinh(image3.image.array[ymin:ymax, xmin:xmax]), origin='lower')
        axs[2].set_title(tag3)
        fig.colorbar(im, ax=axs[2])#, orientation="horizontal")

        outname = "%s/%s-%s-%s_plt.png"%(FIG_FOLDER, tag1, tag2, tag3)
    
    else:
        print("ATT: wrong method!")
        return 1
    
    #plt.tight_layout()
    plt.savefig(outname, bbox_inches="tight", pad_inches=0.1 )

    return 0



def plot_diff_image_catalog(image, catalog):
    fig, ax = plt.subplots(figsize=(7,7))
    display = display_lsp(image, fig, ax)
    display.centroids(catalog, size=20, ctype='r')
    #plt.tight_layout()
    plt.savefig("%s/diff_image_catalog.png"%(FIG_FOLDER), bbox_inches="tight", pad_inches=0.1 )
    
    return 0

    

def plot_crop_grid(image, tag, inj_catalog, crop_size=50, offset=0, if_symmetry=False): 

    len_catalog = len(inj_catalog)
    grid_size = int( len_catalog ** 0.5 )
    big_size = crop_size * grid_size
    big_mat = np.zeros((big_size, big_size))

    big_mat_x_list = []
    big_mat_y_list = []

    for i in range(grid_size):
        for j in range(grid_size):

            ind = i * grid_size + j
            x, y = inj_catalog['x'][ind], inj_catalog['y'][ind]
        
            xmin, xmax = x - crop_size//2, x + crop_size//2
            ymin, ymax = y - crop_size//2, y + crop_size//2

            xmin += offset
            xmax += offset
            ymin += offset
            ymax += offset
        
            big_mat[-(i+1) * crop_size + big_size: -(i) * crop_size + big_size, 
                    j * crop_size: (j+1) * crop_size] = image.image.array[ymin:ymax, xmin:xmax]

            big_mat_x_list.append( (j+0.5) * crop_size )
            big_mat_y_list.append( -(i+0.5) * crop_size + big_size )
            

    plt.figure(figsize=(8,8))
    if if_symmetry:
        threshold = np.max([np.max(big_mat), np.abs(np.min(big_mat))])
        #print("threshold: ", threshold)
        threshold = np.arcsinh(threshold)
        vmin = -threshold
        vmax = threshold
        plt.imshow(np.arcsinh(big_mat), origin="lower", cmap="bwr", vmin=vmin, vmax=vmax)
    else:
        plt.imshow(np.arcsinh(big_mat), origin="lower")

    plt.axis("off")

    plt.plot(big_mat_x_list, big_mat_y_list, "cx")

    plt.savefig("%s/%s_grid.png"%(FIG_FOLDER, tag), bbox_inches="tight", pad_inches=0.1 )
    


#======================================

#def hist_1d(catalog, colname, flag):

#def hist_corner(catalog, colname1, colname2, flag):
