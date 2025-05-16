from lib.tools import *
from lib.imdiff import *


#============================
BAND = 'i'

#============================
diaSources = load_pickle("diaSources")

difference_image = load_pickle("difference_image")

#----------------------------
newline("get_object")
updatedDiaObjects, associatedDiaSources, diaObjectCat, mergedDiaSourceHistory, allDiaSources = get_object(diaSources, difference_image, BAND)
#tmp = get_object(diaSources, difference_image, BAND)
save_pickle("updatedDiaObjects", updatedDiaObjects)
save_pickle("associatedDiaSources", associatedDiaSources)
save_pickle("diaObjectCat", diaObjectCat)
save_pickle("mergedDiaSourceHistory", mergedDiaSourceHistory)
save_pickle("allDiaSources", allDiaSources)

print(diaObjectCat.index.name)
print(diaObjectCat.index.names)

print(allDiaSources.index.name)
print(allDiaSources.index.names)



#----------------------------
newline("get_object2")
# preloadedDiaSources=mergedDiaSourceHistory?

# This gives a RuntimeError: Duplicate DiaSources found after association and merging with history. 
# df.index.has_duplicates
# Note the code is checking index
#updatedDiaObjects, associatedDiaSources, diaObjectCat, mergedDiaSourceHistory, allDiaSources = get_object(diaSources, difference_image, BAND, preloadedDiaSources=allDiaSources, diaObjects=diaObjectCat)

#save_pickle("updatedDiaObjects", updatedDiaObjects)
#save_pickle("associatedDiaSources", associatedDiaSources)
#save_pickle("diaObjectCat", diaObjectCat)