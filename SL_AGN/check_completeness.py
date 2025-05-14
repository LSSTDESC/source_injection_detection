# Load diaSources
# Analyze diaSources
# Associate diaSources

#======================================
from lib.tools import *
from lib.analysis import *

#======================================
inj_catalog = load_pickle("inj_catalog_calexp")
diaSources = load_pickle("diaSources")


#--------------------------------------
plot_tolerance_completeness(inj_catalog, diaSources, diaSrc_tag="base_SdssCentroid")
plot_tolerance_completeness(inj_catalog, diaSources, diaSrc_tag="slot_Centroid")

plot_tolerance_completeness(inj_catalog, diaSources, method="radec")

#plot_completeness_xy(axs[0])
#plot_completeness_xy(axs[1], "base_SdssCentroid")