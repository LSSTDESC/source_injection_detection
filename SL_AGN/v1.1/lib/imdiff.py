#from lib.tools import *

from lsst.ip.diffim.subtractImages import AlardLuptonSubtractTask, AlardLuptonSubtractConfig
from lsst.ip.diffim import detectAndMeasure
#from lsst.ip.diffim import DetectAndMeasureTask #?

#from lsst.ap.association import DiaPipelineTask, DiaPipelineConfig, TransformDiaSourceCatalogTask, AssociationTask
#
#from lsst.meas.base import DetectorVisitIdGeneratorConfig, \
#    DiaObjectCalculationTask

#import pandas as pd



#======================================
# Define some empty tables stored in apdb？
# diaObj
# diaSrc (different from the one generated by diffim)

#DiaObject_empty = pd.read_csv("./DiaObject_empty.csv")
#DiaSource_empty = pd.read_csv("./DiaSource_empty.csv")
#SSObject_empty = pd.read_csv("./SSObject_empty.csv")


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

   

##----------------------------
#def diaSrc2pdDf(diaSourceCat, diffIm, band):
#
#    # Convert LSST AFW SourceCatalog to Pandas dataframe
#    initInputs = {}
#    initInputs['diaSourceSchema'] = diaSourceCat
#    task = TransformDiaSourceCatalogTask(initInputs)
#    res = task.run(diaSourceCat, diffIm, band)
#    
#    return res.diaSourceTable
#    
#    
#
#def get_object(diaSourceCat, diffIm, band, preloadedDiaSources=DiaSource_empty, diaObjects=DiaObject_empty, solarSystemObjectTable=None):
#    
#    # https://pipelines.lsst.io/py-api/lsst.ap.association.DiaPipelineTask.html#lsst.ap.association.DiaPipelineTask.associateDiaSources
#    # https://github.com/lsst/ap_association/blob/7cc9b91cc9da2cce2b05af381bd34bbc5b5e6069/python/lsst/ap/association/diaPipe.py#
#
#    diaSourceTable = diaSrc2pdDf(diaSourceCat, diffIm, band)
#    
#    #diaSourceTable = diaSourceCat.asAstropy().to_pandas()
#    #diaSourceTable['diaSourceId'] = diaSourceTable['id']
#    
#    config = DiaPipelineConfig()
#    config.doConfigureApdb = False
#    #config.apdb.db_url = "sqlite:///file:apdb.db"
#    config.apdb_config_url = "apdb_config.py"
#    
#    dia_pipeline_task = DiaPipelineTask(config=config)
#
##    if diaObjects is None:
##        print("diaObjects is None!")
##        struct = task.createNewDiaObjects(diaSourceTable)
##        print("diaSources: ", struct.diaSources)
##        print("newDiaObjects: ", struct.newDiaObjects)
##        print("nNewDiaObjects: ", struct.nNewDiaObjects)
#
#        #return struct.newDiaObjects
#        
##    else:
#    assocResults = dia_pipeline_task.associateDiaSources(diaSourceTable, solarSystemObjectTable, diffIm, diaObjects)
#    #print(assocResults)
#    #print(len(assocResults))
#    #print(dir(assocResults))
#    #print("associatedDiaSources: ", assocResults.associatedDiaSources)
#    #print("newDiaObjects: ", assocResults.newDiaObjects)
#    #print("newDiaSources: ", assocResults.newDiaSources)
#    #print("marginalDiaSources: ", assocResults.marginalDiaSources)
#        #return struct.newDiaObjects
#
#    associatedDiaSources = assocResults[0] 
#    newDiaObjects = assocResults[1] 
#
#    #------------------------
#    print("\nMerging...")
#    
##    if preloadedDiaSources is None:
##        preloadedDiaSources_new = pd.DataFrame(columns=['diaObjectId', 'diaSourceId'])
##        preloadedDiaSources_new.set_index("diaObjectId", inplace=True)
##        print(preloadedDiaSources_new)
#        
#    mergedDiaSourceHistory, mergedDiaObjects, updatedDiaObjectIds = dia_pipeline_task.mergeAssociatedCatalogs(
#            preloadedDiaSources, associatedDiaSources, 
#            #assocResults.associatedDiaSources,
#            diaObjects, newDiaObjects, 
#            #assocResults.newDiaObjects,
#            diffIm
#        )
#
#    print("\nmergedDiaObjects: ", mergedDiaObjects)
#    
#    #return mergedDiaObjects
#
#    
#    # https://pipelines.lsst.io/py-api/lsst.meas.base.DiaObjectCalculationTask.html#lsst.meas.base.DiaObjectCalculationTask.run
#    dia_object_calculation_task = DiaObjectCalculationTask()
#    diaCalResult = dia_object_calculation_task.run(
#            mergedDiaObjects,
#            mergedDiaSourceHistory,
#            updatedDiaObjectIds,
#            [band])
#
#    print("\ndiaCalResult: ", diaCalResult)
#
#    # Q: diaCalResult.diaObjectCat == mergedDiaObjects? 
#    # No, the latter has not been fully updated
#
#    # preloaded sources should include both {merged history: associated sources + old/preload} and {new src, i.e. unassociated}
#    association_task = AssociationTask()
#    struct = association_task.run(diaSourceTable, diaObjects)
#    matchedDiaSources = struct.matchedDiaSources
#    unAssocDiaSources = struct.unAssocDiaSources
#    nUpdatedDiaObjects = struct.nUpdatedDiaObjects
#    nUnassociatedDiaObjects = struct.nUnassociatedDiaObjects
#
#    
#    allDiaSources = dia_pipeline_task.mergeCatalogs(mergedDiaSourceHistory, unAssocDiaSources, "allDiaSources")
#    
#    allDiaSources.set_index(
#                ["diaObjectId", "band", "diaSourceId"],
#                inplace=True,
#                drop=False)
#    
#    
#    return diaCalResult.updatedDiaObjects, associatedDiaSources, diaCalResult.diaObjectCat, mergedDiaSourceHistory, allDiaSources
