from lib.tools import *
from lib.imdiff import *


#============================
BAND = 'i'

#============================
diaSources = load_pickle("diaSources")

difference_image = load_pickle("difference_image")

#----------------------------
newline("get_object")
updatedDiaObjects, associatedDiaSources, diaObjectCat = get_object(diaSources, difference_image, BAND)
save_pickle("updatedDiaObjects", updatedDiaObjects)
save_pickle("associatedDiaSources", associatedDiaSources)
save_pickle("diaObjectCat", diaObjectCat)

#newDiaObjects = load_pickle("newDiaObjects")
#newDiaObjects2 = get_object(diaSources, difference_image, BAND, diaObjects=newDiaObjects)
#save_pickle("newDiaObjects2", newDiaObjects2)

#----------------------------
