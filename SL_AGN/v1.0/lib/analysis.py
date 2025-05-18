from lib.tools import *

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import hstack
from scipy.spatial import cKDTree

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


#======================================
def match_radec(table_A, table_B, 
                A_ra_col, A_dec_col, A_unit, 
                B_ra_col, B_dec_col, B_unit, 
                tolerance_arcsec=1., if_merge=False): 
    
    # Use Astropy
    # We will use A as the injected catalog and B as the detected catalog

    A_coords = SkyCoord(table_A[A_ra_col], table_A[A_dec_col], unit=A_unit)
    B_coords = SkyCoord(table_B[B_ra_col], table_B[B_dec_col], unit=B_unit)

    # Match A to the nearest neighbor in B
    # Note: it is possible that two A elements matches to one B element.
    # In source injection, if sources in A (injected catalog) are well separated, 
    # then we don't need to worry about it.
    idx, d2d, _ = A_coords.match_to_catalog_sky(B_coords)

    sep_constraint = d2d < (tolerance_arcsec * u.arcsec)

    if if_merge:
        # This can be moved into another function
        A_matches = table_A[sep_constraint]

        B_matches = table_B[idx[sep_constraint]]

        hstack_AB = hstack([A_matches, B_matches])

    else:
        hstack_AB = None
    
    return hstack_AB, sep_constraint, idx[sep_constraint]

    

def match_xy(table_A, table_B, 
             A_x_col, A_y_col, 
             B_x_col, B_y_col, 
             tolerance_pix=1., if_merge=False): 
    
    # Use KDTree
    
    A_coords = np.vstack((table_A[A_x_col], table_A[A_y_col])).T
    B_coords = np.vstack((table_B[B_x_col], table_B[B_y_col])).T

    # Build KDTree for B points
    tree = cKDTree(B_coords)

    # Find nearest neighbor distances for A points
    # distances: An array of minimum distances from each point in A to its nearest neighbor in B.
    # indices: An array of indices in B corresponding to the nearest neighbors.
    # Both distances & indices have the same length as A!
    distances, indices = tree.query(A_coords, k=1)

    matched_A = distances <= tolerance_pix

    if if_merge:
        A_matches = table_A[matched_A]

        B_matches = table_B[indices[matched_A]]

        hstack_AB = hstack([A_matches, B_matches])
    
    else:
        hstack_AB = None
    
    return hstack_AB, matched_A, indices[matched_A]
    

    
#--------------------------------------
def check_completeness(method, 
                       table_A, table_B, 
                       A_col1, A_col2, 
                       B_col1, B_col2, 
                       A_unit=None, B_unit=None, 
                       tolerance=1., 
                       if_update=False, 
                      ):
    
    # Based on the match

    if method=="radec" and A_unit is not None and B_unit is not None:
        A_ra_col, A_dec_col = A_col1, A_col2
        B_ra_col, B_dec_col = B_col1, B_col2

        AB_match, matched_A, B_indices = match_radec(table_A, table_B, 
                                                     A_ra_col, A_dec_col, A_unit, 
                                                     B_ra_col, B_dec_col, B_unit, 
                                                     tolerance_arcsec=tolerance)

    elif method=="xy":
        A_x_col, A_y_col = A_col1, A_col2
        B_x_col, B_y_col = B_col1, B_col2

        AB_match, matched_A, B_indices = match_xy(table_A, table_B, 
                                                  A_x_col, A_y_col, 
                                                  B_x_col, B_y_col, 
                                                  tolerance_pix=tolerance)

    else: 
        print("ATT: wrong method or bad unit!")
        return None, None, None


    cover_rate = np.sum(matched_A) * 1. / len(table_A)
    print("cover_rate: ", cover_rate)

    if if_update:
        B_unique_indices, counts = np.unique(B_indices, return_counts=True)
        B_repeated = np.sum(counts>1)
        print("B_repeated: ", B_repeated )

        RB = np.zeros(len(table_B), dtype=bool)
        RB[B_unique_indices] = True

        # Add RB column to B
        table_B_new = table_B.copy()
        table_B_new["RB"] = RB

    else:
        table_B_new = None

    return cover_rate, table_B_new



#--------------------------------------
def plot_tolerance_completeness(inj_catalog, diaSources, 
                                method="xy", 
                                diaSrc_tag=None, tolerance_pix_max=12):


    # How the completeness evolves with tolerance 
    tolerance_pix = np.arange(0., tolerance_pix_max, 1.)
        
    completeness_list = []

    if method=="xy":
        x_col, y_col = diaSrc_tag + "_x", diaSrc_tag + "_y"
        title = "%s_xy"%diaSrc_tag
    
        for tol in tolerance_pix: 
            print("tol: ", tol)
            completeness = check_completeness(method, 
                                              inj_catalog, diaSources, 
                                              'x', 'y', 
                                              x_col, y_col, 
                                              tol,
                                             )
            
            completeness_list.append(completeness)
            

    elif method=="radec": 
        ra_col, dec_col, unit = "coord_ra", "coord_dec", "rad"
        title = "coord_radec"

        for tol in tolerance_pix: 
            print("tol: ", tol)
            completeness = check_completeness(method, 
                                              inj_catalog, diaSources, 
                                              "ra", "dec", 
                                              ra_col, dec_col, 
                                              "deg", unit, 
                                              tol/ARCSEC2PIX,
                                             )

            completeness_list.append(completeness)
            

    else:
        print("ATT: wrong method!")
        return 1
                                              

    #----------------------------------
    try:
        completeness_is_1_ind = completeness_list.index(1.)
        tolerance_pix_change = tolerance_pix[completeness_is_1_ind]
        print("completeness=1: index, tolerance pix value: ", 
              completeness_is_1_ind, tolerance_pix_change)
    except:
        tolerance_pix_change = None

    
    #----------------------------------
    fig, ax = plt.subplots(1,1)
    
    ax.plot(tolerance_pix, completeness_list, 'o')
    
    ax.axhline(0.5, ls=':', c='grey')
    ax.axhline(1.0, ls=':', c='k')

    if tolerance_pix_change is not None:
        ax.axvline(tolerance_pix[completeness_is_1_ind], ls='--', c='k')

    ax.set_xlabel("Tolerance [pix]")
    ax.set_ylabel("Completeness")
    
    ax.set_title(title)

    ax.set_ylim([-0.1,1.1])
    ax.xaxis.set_major_locator(MultipleLocator(2))  # unit 2
    plt.grid(ls=':', alpha=0.6)

    plt.savefig("%s/tolerance_completeness_%s.png"%(FIG_FOLDER, title), bbox_inches="tight", pad_inches=0.1 )


    return 0


#======================================
    
#======================================        
#def random_forest():