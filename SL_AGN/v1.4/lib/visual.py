import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import lsst.afw.display as afwDisplay

import lib.tools as tl


#============================
afwDisplay.setDefaultBackend("matplotlib")


#============================
def plot_dual(image0, image_inj, x_arr, y_arr, image_type=None):

    if image_type=="visit":
        md_dict = image0.metadata.toDict()
        visit = md_dict["LSST BUTLER DATAID VISIT"]
        visit_tag = "%d"%visit
        tag = "visit_%s"%visit_tag
        
    elif image_type=="template":
        md_dict = image0.metadata.toDict()
        tract = md_dict["LSST BUTLER DATAID TRACT"]
        patch = md_dict["LSST BUTLER DATAID PATCH"]
        template_tag = "%d_%d"%(tract, patch)
        tag = "template_%s"%template_tag
        
    else:
        print("ATT: Wrong image_type!")
        return 1

    
    xlim_min, xlim_max, ylim_min, ylim_max = \
        (
            image0.getBBox().beginX,
            image0.getBBox().endX,
            image0.getBBox().beginY,
            image0.getBBox().endY,
        )

    xmin, xmax = np.min(x_arr), np.max(x_arr)
    ymin, ymax = np.min(y_arr), np.max(y_arr)
    #print("xmin, xmax: ", xmin, xmax)
    #print("ymin, ymax: ", ymin, ymax)

    xlen = xmax - xmin
    ylen = ymax - ymin

    xmin -= 0.05 * xlen
    xmax += 0.05 * xlen
    ymin -= 0.05 * ylen
    ymax += 0.05 * ylen

    if xmin<xlim_min:
        xmin = xlim_min
    if xmax>xlim_max:
        xmax = xlim_max
    if ymin<ylim_min:
        ymin = ylim_min
    if ymax>ylim_max:
        ymax = xlim_max   

    fig, ax = plt.subplots(1, 2, figsize=(10, 5), layout="tight")

    plt.sca(ax[0])
    display1 = afwDisplay.Display(frame=fig)
    display1.scale('asinh', 'zscale')
    display1.image(image_inj.image[xmin:xmax,ymin:ymax])
    plt.title('injected')

    plt.sca(ax[1])
    display2 = afwDisplay.Display(frame=fig)
    display2.scale('asinh', 'zscale')
    display2.image(image0.image[xmin:xmax,ymin:ymax])
    plt.title('original')

    
    plt.savefig('%s/dual_%s.png'%(tl.FIG_FOLDER, tag) )

    return 0