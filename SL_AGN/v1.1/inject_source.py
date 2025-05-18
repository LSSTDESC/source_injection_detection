from lib.tools import *
from lib.inject import *
from lib.butler import *
from lib.visual import *



#============================
SYSTEM_INDEX = 1

TIME_INDEX_LIST = [
0,
]

CALEXP_STAMP_FILENAME_LIST = ["%s/system_%d_%d.fits"%(FIG_FOLDER, SYSTEM_INDEX, i) for i in TIME_INDEX_LIST]
COADD_STAMP_FILENAME = "%s/system_%d_coadd.fits"%(FIG_FOLDER, SYSTEM_INDEX)

MAG_LIST = [
15.,
]

#----------------------------
newline("Start injection!")

print("SYSTEM_INDEX: ", SYSTEM_INDEX)
print("TIME_INDEX_LIST: ", TIME_INDEX_LIST)

print()


#============================
newline("get_calexp_dataId")
calexp_dataId = get_calexp_dataId()

#----------------------------
newline("get_calexp")
calexp = get_calexp(calexp_dataId)
save_pickle("calexp", calexp)

calexp = load_pickle("calexp")
calexp_wcs = calexp.getWcs()

#----------------------------
newline("plot_and_save_image")
tag = "calexp_test"
plot_save_image(calexp, tag)
plot_save_image(calexp, tag, method="lsp")
plot_save_image(calexp, tag, method="plt")

#============================
newline("make_grid_coord")
inj_catalog_calexp = make_grid_coord(calexp_wcs, MAG_LIST[0], CALEXP_STAMP_FILENAME_LIST[0])
print("inj_catalog_calexp: \n", inj_catalog_calexp)

save_pickle("inj_catalog_calexp", inj_catalog_calexp)
inj_catalog_calexp = load_pickle("inj_catalog_calexp")


#============================
newline("calexp_inject_stamp")
injected_calexp = calexp_inject_stamp(calexp, inj_catalog_calexp)
print("injected_calexp: ", injected_calexp)
save_pickle("injected_calexp", injected_calexp)
injected_calexp = load_pickle("injected_calexp")

#----------------------------
plot_save_two_images(calexp, injected_calexp, "calexp", "injected_calexp")



#============================
newline("get_template")
template = get_template(calexp_dataId)
save_pickle("template", template)

template = load_pickle("template")

#----------------------------
inj_catalog_template = inj_catalog_calexp.copy()
inj_catalog_template["stamp"] = [COADD_STAMP_FILENAME] * len(inj_catalog_template)

save_pickle("inj_catalog_template", inj_catalog_template)
inj_catalog_template = load_pickle("inj_catalog_template")


#============================
newline("template_inject_stamp")
injected_template = template_inject_stamp(template, inj_catalog_template)
print("injected_template: ", injected_template)

save_pickle("injected_template", injected_template)
injected_template = load_pickle("injected_template")

#----------------------------
plot_save_two_images(template, injected_template, "template", "injected_template")


#============================
newline("get_src")
sources = get_src(calexp_dataId)

save_pickle("sources", sources)

sources = load_pickle("sources")


#============================
newline("End injection!")
