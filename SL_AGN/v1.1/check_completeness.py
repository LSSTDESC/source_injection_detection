# After DM steps
# Load diaSources
# Analyze diaSources

#======================================
from lib.tools import *
from lib.analysis import *

#======================================
inj_catalog = load_pickle("inj_catalog_calexp")
diaSources = load_pickle("diaSources")
diaSources = diaSources.asAstropy()


#--------------------------------------
plot_tolerance_completeness(inj_catalog, diaSources, diaSrc_tag="base_SdssCentroid")
plot_tolerance_completeness(inj_catalog, diaSources, diaSrc_tag="slot_Centroid")

plot_tolerance_completeness(inj_catalog, diaSources, method="radec")


#--------------------------------------
# Run check_completeness again with a method, a tag, and a fixed tolerance
# The completeness should reach ~1
# Get completeness and marked diaSrc
frac, diaSources_mark = check_completeness(
                            "radec", 
                            inj_catalog, diaSources, 
                            "ra", "dec", 
                            "coord_ra", "coord_dec", 
                            "deg", "rad", 
                            3./ARCSEC2PIX)


# Save diaSrc with R/B marked
# Use load and save pickle
#with open('diaSources_mark.pkl', 'wb') as f:
#    pickle.dump(diaSources_mark, f)
save_pickle("diaSources_mark", diaSources_mark)
diaSources_mark = load_pickle("diaSources_mark")


#--------------------------------------
# Compare columns
# Clearly, bogus have larger error bars on positions!

plot_and_save(diaSources_mark, 'coord_raErr', 'coord_decErr')
plot_and_save(diaSources_mark, 'base_SdssCentroid_xErr', 'base_SdssCentroid_yErr')


# But not covariance
plot_and_save(diaSources_mark, 'coord_ra_dec_Cov', 'coord_ra_dec_Cov')


# The dipolefit columns are not good? Need SNR for flux, and abs() for neg!
plot_and_save(diaSources_mark, 
              'ip_diffim_DipoleFit_pos_instFluxErr', 
              'ip_diffim_DipoleFit_neg_instFluxErr')


#--------------------------------------
# Many other columns good for exploring! Use corner plots.
# May also just look at single columns and their distributions each
# matched and unmatched should have clear separation
# DS tests for example
# Or RF?
