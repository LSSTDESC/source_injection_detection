from lsst.daf.butler import Butler



#============================
BUTLER_CONFIG = "dp02"
COLLECTIONS = "2.2i/runs/DP0.2"
TRACT = 3828
DETECTOR = 19
BAND = 'i'
VISIT_INDEX = 5
BUTLER = Butler(BUTLER_CONFIG, collections=COLLECTIONS) 


#============================
def get_calexp_dataId(butler=BUTLER, tract=TRACT, detector=DETECTOR, band=BAND, visit_index=VISIT_INDEX):

    where = "instrument='LSSTCam-imSim' "
    where += "AND skymap='DC2' "
    where += "AND tract=%d "%tract
    where += "AND detector=%d "%detector
    where += "AND band='%s'"%band
    print("Query: ", where)
    
    calexp_DatasetRefs = sorted(list(set(butler.registry.queryDatasets("calexp", where=where))))
    #print("calexp_DatasetRefs: ", calexp_DatasetRefs)  # very long strings

    # some random visit
    dataId = calexp_DatasetRefs[visit_index].dataId
    #print("calexp_DatasetRefs[%d]: "%VISIT_INDEX, calexp_DatasetRefs[VISIT_INDEX])
    print("dataId: ", dataId)
    print("Type: ", type(dataId) )
    
    return dataId


def get_calexp(calexp_dataId, butler=BUTLER): 
    
    calexp = butler.get("calexp", dataId=calexp_dataId)

    return calexp


def get_src(calexp_dataId, butler=BUTLER): 

    # get sources
    src = butler.get("src", dataId=calexp_dataId)
    
    return src

def get_template(calexp_dataId, butler=BUTLER): 

    # Getting templates/coadds corresponding to the calexp (overlap)
    template = butler.get("goodSeeingDiff_templateExp", dataId=calexp_dataId.required)
    #template = butler.get("deepDiff_templateExp", dataId=calexp_dataId.required)
    
    #template = butler.get("goodSeeingDiff_matchedExp", calexp_dataId) 

    return template





    