from lib.tools import *

from lsst.ip.diffim.subtractImages import AlardLuptonSubtractTask, AlardLuptonSubtractConfig
from lsst.ip.diffim import detectAndMeasure
from lsst.ap.association import DiaPipelineTask

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

   


def get_object(diaSourceTable, diffIm, diaObjects=pd.DataFrame(), solarSystemObjectTable=pd.DataFrame()):
    # https://pipelines.lsst.io/py-api/lsst.ap.association.DiaPipelineTask.html#lsst.ap.association.DiaPipelineTask.associateDiaSources
    
    associatedDiaSources, newDiaObjects = DiaPipelineTask.associateDiaSources(diaSourceTable, solarSystemObjectTable, diffIm, diaObjects)
    print(associatedDiaSources)
    print(newDiaObjects)
    return newDiaObjects
    