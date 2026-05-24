import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import lsst.afw.display as afwDisplay

import lib.tools as tl


#============================
afwDisplay.setDefaultBackend("matplotlib")


#============================
def display_afw(image, fig, ax, title='', scale="asinh"):

    plt.sca(ax)

    display = afwDisplay.Display(frame=fig)
    display.scale(scale, "zscale")
    display.image(image.image)
    
    plt.title(title)

    return display


def plot_two_images(image0, image1, title0, title1, tag):

    fig, axs = plt.subplots(1, 2, figsize=(10, 5), layout="tight")

    display_afw(image0, fig, axs[0], title0)
    display_afw(image1, fig, axs[1], title1)
    
    plt.savefig('%s/%s.png'%(tl.FIG_FOLDER, tag) )
    

def plot_dual(image0, image_inj, inj_radec, tag):
    
    xlim_min, xlim_max, ylim_min, ylim_max = \
        (
            image0.getBBox().beginX,
            image0.getBBox().endX,
            image0.getBBox().beginY,
            image0.getBBox().endY,
        )

    RA_arr = inj_radec['ra']
    DEC_arr = inj_radec['dec']
    wcs = image0.getWcs()
    x_arr, y_arr = wcs.skyToPixelArray(RA_arr, DEC_arr, degrees=True)

    xmin, xmax = np.min(x_arr), np.max(x_arr)
    ymin, ymax = np.min(y_arr), np.max(y_arr)
    #print("xmin, xmax: ", xmin, xmax)
    #print("ymin, ymax: ", ymin, ymax)

    xlen = xmax - xmin
    ylen = ymax - ymin

    xmin -= 0.1 * xlen
    xmax += 0.1 * xlen
    ymin -= 0.1 * ylen
    ymax += 0.1 * ylen

    if xmin<xlim_min:
        xmin = xlim_min
    if xmax>xlim_max:
        xmax = xlim_max
    if ymin<ylim_min:
        ymin = ylim_min
    if ymax>ylim_max:
        ymax = ylim_max

    plot_two_images(image0[xmin:xmax,ymin:ymax], 
                    image_inj[xmin:xmax,ymin:ymax], 
                    "original", "injected", 
                    f"dual_{tag}")

#    fig, ax = plt.subplots(1, 2, figsize=(10, 5), layout="tight")

#    plt.sca(ax[0])
#    display1 = afwDisplay.Display(frame=fig)
#    display1.scale('asinh', 'zscale')
#    display1.image(image_inj.image[xmin:xmax,ymin:ymax])
#    plt.title('injected')

#    plt.sca(ax[1])
#    display2 = afwDisplay.Display(frame=fig)
#    display2.scale('asinh', 'zscale')
#    display2.image(image0.image[xmin:xmax,ymin:ymax])
#    plt.title('original')

#    plt.savefig('%s/dual_%s.png'%(tl.FIG_FOLDER, tag) )

#    return 0


def plot_three_images(image0, image1, image2, title0, title1, title2, tag):

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), layout="tight")

    display_afw(image0, fig, axs[0], title0)
    display_afw(image1, fig, axs[1], title1)
    display_afw(image2, fig, axs[2], title2, "linear")
    
    plt.savefig('%s/%s.png'%(tl.FIG_FOLDER, tag) )



def plot_triple(image0, image1, image2, inj_radec, tag):

    xlim_min, xlim_max, ylim_min, ylim_max = \
        (
            image0.getBBox().beginX,
            image0.getBBox().endX,
            image0.getBBox().beginY,
            image0.getBBox().endY,
        )

    RA_arr = inj_radec['ra']
    DEC_arr = inj_radec['dec']
    wcs = image0.getWcs()
    x_arr, y_arr = wcs.skyToPixelArray(RA_arr, DEC_arr, degrees=True)

    xmin, xmax = np.min(x_arr), np.max(x_arr)
    ymin, ymax = np.min(y_arr), np.max(y_arr)

    xlen = xmax - xmin
    ylen = ymax - ymin

    xmin -= 0.1 * xlen
    xmax += 0.1 * xlen
    ymin -= 0.1 * ylen
    ymax += 0.1 * ylen

    if xmin<xlim_min:
        xmin = xlim_min
    if xmax>xlim_max:
        xmax = xlim_max
    if ymin<ylim_min:
        ymin = ylim_min
    if ymax>ylim_max:
        ymax = ylim_max   

    plot_three_images(image0[xmin:xmax,ymin:ymax], 
                      image1[xmin:xmax,ymin:ymax], 
                      image2[xmin:xmax,ymin:ymax], 
                      "visit", "warped_template", "difference",  
                      f"triple_{tag}")