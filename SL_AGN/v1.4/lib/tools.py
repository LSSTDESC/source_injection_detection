import numpy as np
import pickle


#======================================
FIG_FOLDER = "fig"
STAMP_FOLDER = "stamp"

CATALOG_FOLDER = "catalog"

ARCSEC2PIX = 0.2


#======================================
def mag2flux(mag, m0):
    return 10**( (mag - m0)/(-2.5) )


def flux2mag(flux, m0):
    return -2.5 * np.log10(flux) + m0


#======================================
def newline(prompt=""):
    print('\n' + '-'*20)
    print(prompt)
    print('-'*20)


def save_pickle(tag, item, folder):
    filename = "%s/%s.pkl"%(folder, tag) 
    with open(filename, "wb") as f:
        pickle.dump(item, f)

    return 0
    

def load_pickle(tag):
    filename = "%s/%s.pkl"%(folder, tag) 
    with open(filename, "rb") as f:
        item = pickle.load(f)

    return item


#======================================
def get_visit_rotation_angle(visit): 
    
    rotation_angle0 = visit.visitInfo.getBoresightRotAngle()
    rotation_angle = -rotation_angle0.asDegrees()
    
    return rotation_angle