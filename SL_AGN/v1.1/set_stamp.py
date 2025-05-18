from lib.tools import *
from lib.stamp import *



#============================
SYSTEM_INDEX = 1

TIME_INDEX_LIST = [
0,
]


#----------------------------
newline("Start!")

print("SYSTEM_INDEX: ", SYSTEM_INDEX)
print("TIME_INDEX_LIST: ", TIME_INDEX_LIST)

print()


#============================
newline("get_single_stamp")
get_single_stamp(SYSTEM_INDEX, TIME_INDEX_LIST[0])

#----------------------------
newline("get_coadd_stamp")
get_coadd_stamp(SYSTEM_INDEX)

#----------------------------
newline("get_diff_stamp")
get_diff_stamp(SYSTEM_INDEX, TIME_INDEX_LIST[0])

#----------------------------
newline("check_flux_diff")
check_flux_diff(SYSTEM_INDEX)

#----------------------------
newline("add_wcs")
filename = "%s/system_%d_%d.fits"%(FIG_FOLDER, SYSTEM_INDEX, TIME_INDEX_LIST[0])
add_wcs(filename)

filename = "%s/system_%d_coadd.fits"%(FIG_FOLDER, SYSTEM_INDEX)
add_wcs(filename)


#============================
newline("End!")
