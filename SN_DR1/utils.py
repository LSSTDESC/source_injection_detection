import numpy as np
from astropy.time import Time
from astropy.table import Table
import pandas as pd
import random
import sncosmo
from scipy.special import gamma

# LSST-specific things
from lsst.rsp import get_tap_service
from lsst.daf.butler import Butler, Timespan
import lsst.geom as geom
from lsst.source.injection import VisitInjectConfig, VisitInjectTask

# Import the services
service = get_tap_service("tap")
assert service is not None

butler = Butler("dp1", collections="LSSTComCam/DP1")
assert butler is not None

pixelscale = 0.2
bands = ["u", "g", "r", "i", "z", "y"]

def sersic_n_to_b(n):
    """Compute the `b(n)` for a sersic model. This factor ensures that
    the $R_e$ and $I_e$ parameters do in fact correspond
    to the half light values and not some other scale
    radius/intensity.

    """

    return (
        2 * n
        - 1 / 3
        + 4 / (405 * n)
        + 46 / (25515 * n**2)
        + 131 / (1148175 * n**3)
        - 2194697 / (30690717750 * n**4)
    )

def sersic_flux_to_Ie(flux, n, Re, q):
    """Compute the intensity at $R_e$ (flux/arcsec^2) for a 2D
    elliptical sersic given the $F,n,R_e,q$ parameters which
    uniquely define the profile ($F$ is the total flux
    integrated to infinity). Note that $R_e$ is the effective
    radius in the sersic representation:

    $$I(R) = I_ee^{-b_n[(R/R_e)^{1/n}-1]}$$

    **Args:**
    -  `flux`: flux integrated to infinity (flux)
    -  `n`: sersic index
    -  `Re`: Radius at half integrated light
    -  `q`: axis ratio (b/a)

    """
    bn = sersic_n_to_b(n)
    return flux / (2 * np.pi * Re**2 * q * n * (np.exp(bn) * bn ** (-2 * n)) * gamma(2 * n))

def gather_host_data(objectID, saveto="host_data.npz"):
    query = (
        "SELECT coord_ra, coord_dec, ObjectId, refExtendedness, refSizeExtendedness, sersic_ra, sersic_dec, sersic_reff_x, sersic_reff_y, sersic_index, sersic_rho, u_sersicFlux, g_sersicFlux, r_sersicFlux, i_sersicFlux, z_sersicFlux, y_sersicFlux "
        "FROM dp1.Object "
        f"WHERE ObjectId = {objectID}"
    )
    
    job = service.submit_job(query)
    job.run()
    job.wait(phases=["COMPLETED", "ERROR"])
    print("Job phase is", job.phase)
    job.raise_if_error()
    assert job.phase == "COMPLETED"
    
    df_host = job.fetch_result().to_table()
    print(f"\nFound {len(df_host)} source(s)")
    assert len(df_host) == 1

    host_data = {
        "u": {},
        "g": {},
        "r": {},
        "i": {},
        "z": {},
        "y": {},
        "ra": float(df_host["sersic_ra"][0]),
        "dec": float(df_host["sersic_dec"][0])
    }
    for band in bands:
        dr2 = df_host["sersic_reff_x"][0]**2 - df_host["sersic_reff_y"][0]**2
        ac = 2 * df_host["sersic_rho"][0] * df_host["sersic_reff_x"][0] * df_host["sersic_reff_y"][0]
        q_1 = df_host["sersic_reff_x"][0]**2 + df_host["sersic_reff_y"][0]**2
        q_2 = np.sqrt(dr2**2 + ac**2)
        host_data[band]["q"]   = float(np.sqrt((q_1 - q_2) / (q_1 + q_2)))
        host_data[band]["PA"]  = float(0.5 * np.arctan2(ac, dr2))
        host_data[band]["n"]   = float(df_host["sersic_index"][0])
        host_data[band]["Re"]  = float(np.sqrt(q_1) * pixelscale)
        host_data[band]["Ie"]  = float(sersic_flux_to_Ie(df_host[f"{band}_sersicFlux"][0], host_data[band]["n"], host_data[band]["Re"], host_data[band]["q"]))
    np.savez(saveto, **host_data)

##### Code useful for Part 2 of the plan #####
def dataset_ref_to_dataframe(dataset_refs):
    """
    Convert DatasetRef object(s) to a pandas DataFrame.
    
    Parameters
    ----------
    dataset_refs : DatasetRef or list of DatasetRef
        A single Butler DatasetRef object or a list of them.
    
    Returns
    -------
    pd.DataFrame
        A DataFrame containing all the dataset reference information.
        Each row represents one DatasetRef.
    """
    # Handle single DatasetRef by converting to list
    if not isinstance(dataset_refs, list):
        dataset_refs = [dataset_refs]
    
    # Collect data from all DatasetRefs
    data = []
    for dataset_ref in dataset_refs:
        # Extract dataset type information
        dataset_type = dataset_ref.datasetType
        
        # Convert dataId to a dictionary
        data_id_dict = dict(dataset_ref.dataId.mapping)
        
        # Build a dictionary with all the information
        info_dict = {**data_id_dict,}
        
        data.append(info_dict)
    
    return pd.DataFrame(data)

def find_visits_from_ObjectId(
    objId:int=611255003922842442,
    min_time:float=45_000.,# Corresponds to year ~1982.
    max_time:float=69_807.,# Correspond to year ~2050.
    timespan=None,
    ):
    """
    Finds the visits (and some parameters of the visits, like the filter, day of observation,
    etc.) from a given ObjectId, a minimum observation time, and a maximum observation time.

    Input:
    - objId: int: Object ID.
    - min_time: float: Minimum time of observation in MJD.
    - max_time: float: Maximum time of observation in MJD.

    Returns:
    - dictionary of all the visit IDs and visit times.
    """
    # times of interest
    if timespan is None:
        time1 = Time(min_time, format="mjd", scale="tai")
        time2 = Time(max_time, format="mjd", scale="tai")  
        timespan = Timespan(time1, time2)

    # First, start with querying the Object database for the object's coordinates:
    job = service.submit_job(
        "SELECT obj.objectId, obj.coord_ra, obj.coord_dec " + \
        "FROM dp1.Object AS obj " + \
        f"WHERE obj.objectId = {objId}"
    )
    job.run()
    job.wait(phases=['COMPLETED', 'ERROR'])
    job.raise_if_error()
    assert job.phase == 'COMPLETED'
    
    # Fetch the results (ra and dec of object):
    results = job.fetch_result()
    tab = results.to_table()
    ra, dec = tab['coord_ra'][0], tab['coord_dec'][0]
    del results, tab

    # Now, we make a new query on the visits database:
    query = f"visit_detector_region.region OVERLAPS POINT({ra}, {dec}) " + \
            "AND band.name IN ('u', 'g', 'r', 'i', 'z', 'y') " + \
            "AND visit.timespan OVERLAPS :timespan "
    visit_img_refs = butler.query_datasets(
        'visit_image', where=query, bind={"timespan": timespan}, order_by='visit.timespan.begin', with_dimension_records=True
    )
    # Transform the reference dataset into a pandas DataFrame:
    visit_img_df = dataset_ref_to_dataframe(visit_img_refs)

    # Get exact visit times
    job = service.submit_job(
        "SELECT vis.visit, vis.expMidptMJD "
        "FROM dp1.Visit as vis "
        f"WHERE vis.visit IN ({', '.join(list(str(v) for v in visit_img_df['visit']))})"
    )
    job.run()
    job.wait(phases=["COMPLETED", "ERROR"])
    job.raise_if_error()
    assert job.phase == "COMPLETED"
    df_exposure = job.fetch_result().to_table().to_pandas()

    # Add MJD column
    visit_img_df = visit_img_df.merge(df_exposure, on="visit", how="left")
    return visit_img_df, ra, dec    

##### Code useful for Part 3 of the plan #####

def asymgaus(mean, sigma_m, sigma_p):
    """
    Returns an asymmetric Gaussian distribution.
    """
    return mean + np.random.choice([sigma_m, sigma_p])*np.random.normal()

def generate_SN(z, t0, band='lssti', source='salt3',):
    """
    Returns an sncosmo.Model object containing all the supernova intrinsic parameters:
    - x0: derived from intrinsic absolute magnitude.
    - c, x1: from Scolnic & Kessler.
    - z: to be set by the host galaxy.
    - t0: to be set by the higher-order function which calls generate_SN.
    """
    while True:
        x1 = asymgaus(0.973, 1.472, 0.222) # Scholnic and Kessler, Table 1, Var: G10, Survey: high z
        if -3 < x1 < 2:
            break
    while True:
        c = asymgaus(-0.054, 0.043, 0.101) # Scholnic and Kessler, Table 1, Var: G10, Survey: high z
        if -0.3 < c < 0.5:
            break
    while True:
        M = np.random.normal(loc=-19.35, scale=0.13)  # Gaussian scatter on M
        if -20 < M < -19:
            break
    sig_coh = np.random.normal(scale=0.09)
    sig_incoh = np.random.normal(scale=0.07)
    M = M + sig_coh + sig_incoh
    sn = sncosmo.Model(source=source)
    sn.set( # Set the SN parameters EXCEPT the absolute magnitude:
        z = z,
        t0 = t0,
        x1 = x1,
        c = c,
        )
    sn.set_source_peakabsmag(M, band, 'ab') # Set absolute mangitude:
    return sn

def get_an_SN_and_object(
    visit_info,
    z=0.1,# Redshift of the galaxy (used to set the SNIa redshift)
    t0=60_635.246,# MJD value for the intrinsic max flux of the SNIa in the i-band.
    ):
    """
    
    """
    visit_info = visit_info.copy()
    
    # Save the number of observations:
    num_obs = len(visit_info)
    if num_obs == 0:
        print("There are no observations for this object")
        return None
    
    # Then, generate the desired SNIa:
    sn = generate_SN(z=z, t0=t0)

    # Create the observation table as an astropt.table.Table:
    obs_table = Table({
        'time': visit_info['expMidptMJD'].values,
        'band': ['lsst' + k for k in visit_info['band'].values],
        'gain': [1.] * num_obs,
        'skynoise': [0.] * num_obs,
        'zp': [31.4] * num_obs,
        'zpsys':['ab'] * num_obs,
    })
    
    # Instantiate the model parameters and create the light curve:
    model_params = {name: sn.get(name) for name in sn.param_names}
    lcs = sncosmo.realize_lcs(obs_table, sn, params=[model_params], scatter=False)

    # Add flux array as new column in the galaxy_observations DataFrame:
    visit_info['flux'] = lcs[0]['flux']
    return visit_info

def flux_to_mag(flux, band):
    mag = -2.5 * np.log10(flux) + 31.4
    return mag

def inject_visit_images(visit_info, ra, dec, cutout_size=128):
    injected_data = []
    
    num_obs = len(visit_info)
    for i in range(num_obs):
        # get visit image
        dataId = {
            "instrument": visit_info.at[i, "instrument"], 
            "detector":  visit_info.at[i, "detector"], 
            "visit":  visit_info.at[i, "visit"], 
            "band":  visit_info.at[i, "band"], 
            "day_obs":  visit_info.at[i, "day_obs"], 
            "physical_filter":  visit_info.at[i, "physical_filter"]
        }
        visit_img = butler.get('visit_image', dataId=dataId)

        # Build injection info
        my_injection_catalog = Table(
            {
                'injection_id': [9999],
                'ra': [ra],
                'dec': [dec],
                "seed": ["123"],
                'source_type': ['Star'],
                'mag': [flux_to_mag(visit_info.at[i, "flux"], visit_info.at[i, "band"])],
            }
        )

        # Run injection
        inject_config = VisitInjectConfig()
        inject_task = VisitInjectTask(config=inject_config)
        injected_output = inject_task.run(
            injection_catalogs=my_injection_catalog,
            input_exposure=visit_img.clone(),
            psf=visit_img.psf.clone(),
            photo_calib=visit_img.photoCalib,
            wcs=visit_img.wcs,
        )

        # Get cutout
        try:
            injected_exposure = injected_output.output_exposure
            spherePoint = geom.SpherePoint(ra*geom.degrees, dec*geom.degrees)
            cutout = injected_exposure.getCutout(spherePoint, geom.Extent2I(cutout_size))
            mask_list = ["BAD", "CLIPPED", "CR", "CROSSTALK", "EDGE", "INTRP", "NO_DATA", "REJECTED", "SAT", "SENSOR_EDGE", "STREAK", "SUSPECT", "UNMASKEDNAN", "VIGNETTED"]
            combined_bitmask = 0
            for plane in mask_list:
                combined_bitmask |= cutout.mask.getPlaneBitMask(plane)
            psf_obj = cutout.getPsf() # fixme maybe psf is more stable from visit image?
            xy = cutout.wcs.skyToPixelArray(ra, dec, degrees=True)
            point = geom.Point2D(xy[0].item(), xy[1].item())
            psf_img = psf_obj.computeImage(point)
            offset = geom.Extent2D(geom.Point2I(0, 0) - cutout.getXY0())
            shifted_wcs = cutout.getWcs().getFitsApproximation().copyAtShiftedPixelOrigin(offset)
        except Exception as e:
            print("Encountered ERROR: ", e)
            injected_data.append(None)
            continue

        # save results
        injected_data.append({
            "image": cutout.image.array,
            "flux": visit_info.at[i, "flux"],
            "mjd": visit_info.at[i, 'expMidptMJD'],
            "variance": cutout.variance.array,
            "mask": (cutout.mask.array & combined_bitmask) != 0,
            "psf": psf_img.array,
            "fits metadata": dict(shifted_wcs.getFitsMetadata()),
            **dataId,
        })
    return injected_data

def injected_data_to_sn_data(injected_data, ra, dec, saveto = "sn_data.npz"):
    sn_data = {
        "ra": ra,
        "dec": dec,
        "image": [],
        "flux": [],
        "mjd": [],
        "variance": [],
        "mask": [],
        "psf": [],
        "fits metadata": [],
        "instrument": [],
        "detector": [],
        "visit": [],
        "band": [],
        "day_obs": [],
        "physical_filter": [],
    }

    injected_data = list(filter(lambda d: d is not None, injected_data))
    mjd = np.array(list(inj_dat["mjd"] for inj_dat in injected_data))
    N = np.argsort(mjd)
    for i in N:
        dat = injected_data[i]
        for key in dat:
            sn_data[key].append(dat[key])
    for key in sn_data:
        if key in ["ra", "dec", "fits metadata"]:
            continue
        sn_data[key] = np.stack(sn_data[key])
    np.savez(saveto, **sn_data)

if __name__=='main':
    find_visits_from_ObjectId()