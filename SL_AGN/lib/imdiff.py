from lib.tools import *

from lsst.ip.diffim.subtractImages import AlardLuptonSubtractTask, AlardLuptonSubtractConfig
from lsst.ip.diffim import detectAndMeasure
#from lsst.ip.diffim import DetectAndMeasureTask #?

from lsst.ap.association import DiaPipelineTask, DiaPipelineConfig, TransformDiaSourceCatalogTask

import pandas as pd


#======================================
def run_AL(templateExposure, scienceExposure, sources): 

    # image subtraction
    
    config = AlardLuptonSubtractConfig()
    
    try:
        config.sourceSelector.value.unresolved.name = "base_ClassificationExtendedness_value" # 1 for extended source
    except:
        pass

    alTask = AlardLuptonSubtractTask(config=config)

    al_result = alTask.run(templateExposure, scienceExposure, sources)

    return al_result



def detect_measure(scienceExposure, templateExposure, difference): 
    
    # running detection on a difference image (al_result.difference)

    task = detectAndMeasure.DetectAndMeasureTask()
    
    diff_dm_result = task.run(scienceExposure, templateExposure, difference)

    return diff_dm_result

   

def diaSrc2pdDf(diaSourceCat, diffIm, band):

    # Convert SourceCatalog to Pandas dataframe
    initInputs = {}
    initInputs['diaSourceSchema'] = diaSourceCat
    task = TransformDiaSourceCatalogTask(initInputs)
    res = task.run(diaSourceCat, diffIm, band)
    
    return res.diaSourceTable
    
    

def get_object(diaSourceCat, diffIm, band, diaObjects=None, solarSystemObjectTable=None):
    
    # https://pipelines.lsst.io/py-api/lsst.ap.association.DiaPipelineTask.html#lsst.ap.association.DiaPipelineTask.associateDiaSources
    # https://github.com/lsst/ap_association/blob/7cc9b91cc9da2cce2b05af381bd34bbc5b5e6069/python/lsst/ap/association/diaPipe.py#

    # It looks like we need to merge the new DIA object catalog with the old one
    # Maybe it's better to run the dia pipeline as a whole

    diaSourceTable = diaSrc2pdDf(diaSourceCat, diffIm, band)
    
    config = DiaPipelineConfig()
    config.apdb_config_url = "apdb_config.py"
    
    task = DiaPipelineTask(config=config)

    if diaObjects is None:
        print("diaObjects is None!")
        struct = task.createNewDiaObjects(diaSourceTable)
        print("diaSources: ", struct.diaSources)
        print("newDiaObjects: ", struct.newDiaObjects)
        print("nNewDiaObjects: ", struct.nNewDiaObjects)

        return struct.newDiaObjects
        
    else:
        struct = task.associateDiaSources(diaSourceTable, solarSystemObjectTable, diffIm, diaObjects)
        print("associatedDiaSources: ", struct.associatedDiaSources)
        print("newDiaObjects: ", struct.newDiaObjects)
        print("newDiaSources: ", struct.newDiaSources)
        print("marginalDiaSources: ", struct.marginalDiaSources)
        return struct.newDiaObjects
    