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
from matplotlib.ticker import ScalarFormatter

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

import pandas as pd


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
        return None, None


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
            #print("tol: ", tol)
            completeness, _ = check_completeness(method, 
                                                 inj_catalog, diaSources, 
                                                 'x', 'y', 
                                                 x_col, y_col, 
                                                 tol,
                                                 )
            #print("completeness: ", completeness)
            
            completeness_list.append(completeness)
            

    elif method=="radec": 
        ra_col, dec_col, unit = "coord_ra", "coord_dec", "rad"
        title = "coord_radec"

        for tol in tolerance_pix: 
            #print("tol: ", tol)
            completeness, _ = check_completeness(method, 
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
    fig, ax = plt.subplots(1,1,layout="constrained",figsize=(4,3))

    #print("tolerance_pix, completeness_list: \n", tolerance_pix, completeness_list)
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
#def remove_nans(data):
#    return data[~np.isnan(data).any(axis=1)]


def plot_corner(B, x_col, y_col, bins=20, save_dir=FIG_FOLDER):

    matched = B[B["RB"] == True]
    unmatched = B[B["RB"] == False]

    x_min = min(np.nanmin(matched[x_col]), np.nanmin(unmatched[x_col]))
    x_max = max(np.nanmax(matched[x_col]), np.nanmax(unmatched[x_col]))
    x_bin_edges = np.linspace(x_min, x_max, bins + 1)
    x_inv = x_max - x_min
    x_min -= x_inv * 0.1
    x_max += x_inv * 0.1
    if np.isnan(x_min) or np.isnan(x_max):
        return 1
    
    y_min = min(np.nanmin(matched[y_col]), np.nanmin(unmatched[y_col]))
    y_max = max(np.nanmax(matched[y_col]), np.nanmax(unmatched[y_col]))
    y_bin_edges = np.linspace(y_min, y_max, bins + 1)
    y_inv = y_max - y_min
    y_min -= y_inv * 0.1
    y_max += y_inv * 0.1
    if np.isnan(y_min) or np.isnan(y_max):
        return 1
    
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(4, 4),
                           gridspec_kw={'width_ratios': [1, 0.2], 
                                    'height_ratios': [0.2, 1]},
                           layout="constrained",
                           #sharex=True, sharey=True
                          )

    
# Scatter plot in the main axes (bottom-left)
    ax[1, 0].scatter(matched[x_col], matched[y_col], color='b', label='True', alpha=0.2, marker='^')
    ax[1, 0].scatter(unmatched[x_col], unmatched[y_col], color='r', label='False', alpha=0.2)
    ax[1, 0].set_xlabel(x_col)
    ax[1, 0].set_ylabel(y_col)
    ax[1, 0].set_xlim(x_min, x_max)
    ax[1, 0].set_ylim(y_min, y_max)
    ax[1, 0].legend()

# 2D histogram in the main axes (bottom-left)
    #hb = ax[1, 0].hist2d(data['x'], data['y'], bins=30, cmap='Blues', alpha=0.5)
    #fig.colorbar(hb[3], ax=ax[1, 0])

# Histograms of x on top (top-left)
    ax[0, 0].hist(matched[x_col], bins=x_bin_edges, color='b', alpha=0.7, label='True')#, log=True)
    ax[0, 0].hist(unmatched[x_col], bins=x_bin_edges, color='r', alpha=0.7, label='False')#, log=True)
    #ax[0, 0].set_ylabel('Frequency')
    #ax[0, 0].legend()
    ax[0, 0].set_xlim(x_min, x_max)
    ax[0, 0].set_xticks([])
    ax[0, 0].axvline(np.median(matched[x_col]), ls=':', c='b')
    ax[0, 0].axvline(np.median(unmatched[x_col]), ls=':', c='r')
    print("line: ", np.median(matched[x_col]), np.median(unmatched[x_col]), )
    

# Histograms of y on the right (right-top)
    ax[1, 1].hist(matched[y_col], bins=y_bin_edges, color='b', alpha=0.7, label='True', orientation='horizontal')#, log=True, )
    ax[1, 1].hist(unmatched[y_col], bins=y_bin_edges, color='r', alpha=0.7, label='False', orientation='horizontal')#, log=True, )
    #ax[1, 1].set_xlabel('Frequency')
    #ax[1, 1].legend()
    ax[1, 1].set_ylim(y_min, y_max)
    ax[1, 1].set_yticks([])
    ax[1, 1].axhline(np.median(matched[y_col]), ls=':', c='b')
    ax[1, 1].axhline(np.median(unmatched[y_col]), ls=':', c='r')
    print("line: ", np.median(matched[y_col]), np.median(unmatched[y_col]), )

    ax[0, 1].axis('off')
# Adjust layout to make room for the top and right histograms
    #plt.tight_layout()
    #plt.show()

    formatter = ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, 3))  # Control the threshold for using scientific notation

    # Apply the formatter to both x and y axes
    ax[1, 0].xaxis.set_major_formatter(formatter)
    ax[1, 0].yaxis.set_major_formatter(formatter)

    #print(ax[1, 0].get_yticklabels())

    #ax[1, 0].yaxis.offsetText.set_visible(False)
    #print(ax[1, 0].yaxis.offsetText.get_text())
    #offset = ax[1, 0].yaxis.get_major_formatter().get_offset()

    #print(ax[1, 0].yaxis.offsetText)
    #print(dir(ax[1, 0]))
    #ax[1, 0].yaxis.set_label_text(ax[1, 0].get_ylabel() + " " + offset)

    filename = f"{save_dir}/{x_col}_{y_col}.png"

    # Adjust layout and save figure
    #plt.tight_layout()
    #fig.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.savefig(filename) #, dpi=300)
    print("savefig: ", filename)
    

    return 0


    
#======================================        
def random_forest(diaSources_mark, col_list):
    

    diaSources_mark_s = diaSources_mark[col_list]
    diaSources_mark_s_arr = np.array(diaSources_mark_s.to_pandas().values, dtype=np.float32)

    # X and y for input values and flags
    X = diaSources_mark_s_arr[:,:-1]
    y = diaSources_mark_s_arr[:,-1]

    print("X[:3]: \n", X[:3])
    print("y[:3]: \n", y[:3])

    #------------------------
    # Divide the sample into training and test datasets
    X_train, X_test, y_train, y_test = \
                    train_test_split(X, y, 
                                     test_size=0.2, 
                                     random_state=42)

    #------------------------
    # Training
    rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_classifier.fit(X_train, y_train)

    #------------------------
    #Predict
    y_pred = rf_classifier.predict(X_test)

    accuracy = accuracy_score(y_test, y_pred)

    print(f"Accuracy: {accuracy:.2f}")

    return 0
