from lib.tools import *
from lib.imdiff import *


#======================================
diaSources = load_pickle("diaSources")
difference_image = load_pickle("difference_image")
newDiaObjects = get_object(diaSource, difference_image)
