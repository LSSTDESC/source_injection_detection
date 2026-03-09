from lsst.daf.butler import Butler



#============================
BUTLER_CONFIG = "dp1"
COLLECTIONS = "LSSTComCam/DP1"
BUTLER = Butler(BUTLER_CONFIG, collections=COLLECTIONS) 


#============================
def get_visit_dataset_refs(ra, dec, band, butler=BUTLER):

    dataset_refs = butler.query_datasets("visit_image",
                                     where="band.name = :band AND \
                                     visit_detector_region.region OVERLAPS POINT(:ra, :dec)",
                                     bind={"band": band,
                                           "ra": ra, "dec": dec},
                                     order_by=["visit.timespan.begin"])
    
    return dataset_refs


def get_visit_image(dataset_refs, visit_index, butler=BUTLER): 
    
    ref = dataset_refs[visit_index]
    visit_image = butler.get(ref)

    return visit_image


def get_template_dataset_refs(ra, dec, band, butler=BUTLER): 

    query = f"band.name = '{band}' AND patch.region OVERLAPS POINT({ra}, {dec})"
    dataset_refs = butler.query_datasets("template_coadd", where=query)

    return dataset_refs


def get_template_image(dataset_refs, template_index, butler=BUTLER):
    
    ref = dataset_refs[template_index]
    template_image = butler.get(ref)

    return template_image




    