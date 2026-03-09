# generate_slsim_lsst_dataset.py

### Cosmology and astropy packages
from astropy.cosmology import FlatLambdaCDM
from astropy.units import Quantity
import astropy.coordinates as coord
import astropy.units as u

### Arrays, tables, plots
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy

### SLSim functions
import slsim.Sources as sources
import slsim.Deflectors as deflectors
import slsim.Pipelines as pipelines
from slsim.Sources.SourceCatalogues.QuasarCatalog.quasar_pop import QuasarRate
from slsim.Lenses.lens_pop import LensPop
from slsim.ImageSimulation.image_simulation import (
    point_source_coordinate_properties,
    lens_image_series_precomputed_mags,
    lens_image
)
import h5py

from slsim.Util.param_util import ellipticity2phi_q

# from slsim.Util.distribution_plot_utils import make_contour
from slsim.LsstSciencePipeline.rubin_sim_pipeline import get_rubin_cadence
from lenstronomy.Util.data_util import bkg_noise

### Readin, readout, paths
from contextlib import redirect_stdout
import io
from tqdm import tqdm
import argparse
import json
from datetime import datetime

lsst_colors = {
    "u": "#0c71ff",
    "g": "#49be61",
    "r": "#c61c00",
    "i": "#ffc200",
    "z": "#f341a2",
    "y": "#5d0000",
}

def get_random_ra_dec(N=1):
        ra_points = coord.Angle(np.random.uniform(low=0, high=360, size=N) * u.degree)
        ra_points = ra_points.wrap_at(180 * u.degree)
        # dec goes from -72 to +12
        lower = -70
        upper = 10
        p = (
            np.sin(np.random.uniform(low=lower, high=upper, size=N) * u.deg)
            - np.sin(lower * u.deg)
        ) / (np.sin(upper * u.deg) - np.sin(lower * u.deg))
        dec_points = coord.Angle(
            ((((np.arcsin(2 * p - 1).to(u.deg) + 90 * u.deg) / (180 * u.deg)) * 84) + lower)
            * u.deg
        )
        return ra_points, dec_points

def compute_magnitude_zeropoint(mag_zp_1s, exposure_time=30, gain=1):
        return mag_zp_1s + 2.5 * np.log10(exposure_time / gain)

class DatasetGenerator:
    def __init__(self, sky_area, bands, postage_stamps_path, dataset_path):
        self.bands = bands
        self.full_sky_area = sky_area
        self.postage_stamps_path = postage_stamps_path
        self.dataset_path = dataset_path


    def set_up_cosmology(self):
        # define a cosmology
        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

        # define a sky area
        self.galaxy_sky_area = Quantity(
            value=10, unit="deg2"
        )  # this is the sky area over which galaxies are sampled
        self.quasar_sky_area = Quantity(value=10, unit="deg2")

        # this is the sky area over which lensed quasars are sampled
        self.full_sky_area = Quantity(value=self.full_sky_area, unit="deg2")

        # define limits in the intrinsic deflector and source population (in addition
        # to the skypy config
        # file)
        # if this is to faint, you waste computing
        
        self.kwargs_deflector_cut = {"band": "i", "band_max": 28, "z_min": 0.01, "z_max": 2.5}
        self.kwargs_source_cut = {"band": "i", "band_max": 26, "z_min": 0.001, "z_max": 6.0}

    # generate galaxy population using skypy pipeline.
    def set_up_galaxy_quasar_population(self):
        self.galaxy_simulation_pipeline = pipelines.SkyPyPipeline(
        skypy_config=None,
        sky_area=self.galaxy_sky_area,
        filters=self.bands,
        cosmo=self.cosmo,
        z_min=0,
        )
        
        self.lens_galaxies_ell = deflectors.EllipticalLensGalaxies(
        galaxy_list=self.galaxy_simulation_pipeline.red_galaxies,
        kwargs_cut=self.kwargs_deflector_cut,
        kwargs_mass2light={},
        cosmo=self.cosmo,
        sky_area=self.galaxy_sky_area,
        gamma_pl=dict(mean=2.0, std_dev=0.16),
        )

        # Initiate QuasarRate class to generate quasar sample.
        self.quasar_class = QuasarRate(
            cosmo=self.cosmo,
            sky_area=self.quasar_sky_area,
            noise=True,
            redshifts=np.linspace(0.001, 6.00, 100),  # these redshifts are provided
            # to match general slsim redshift range in skypy pipeline.
        )
        # quasar sample with host galaxy
        self.quasar_source_plus_galaxy = self.quasar_class.quasar_sample(
            m_min=15, m_max=28, host_galaxy=True
        )

        # Prepare dictionary of agn variability kwargs
        length_of_light_curve = 3850

        # log(BH_mass/Msun), M_i, log(SFi_inf/mag), log(tau/days), zsrc
        MACLEOD2010_MEANS = np.array(
            [8.53308079, -23.48721021, -0.51665998, 2.28708691, 2.11640976]
        )
        MACLEOD2010_COV = np.array(
            [
                [0.27862905, -0.29501766, 0.00675703, 0.04606804, -0.00665875],
                [-0.29501766, 2.06855169, 0.19690851, 0.0244139, -0.29913764],
                [0.00675703, 0.19690851, 0.02785685, 0.01083628, -0.02216221],
                [0.04606804, 0.0244139, 0.01083628, 0.05636087, -0.02716507],
                [-0.00665875, -0.29913764, -0.02216221, -0.02716507, 0.3077278],
            ]
        )

        variable_agn_kwarg_dict = {
            "multivariate_gaussian_means": MACLEOD2010_MEANS,
            "multivariate_gaussian_covs": MACLEOD2010_COV,
            "known_band": "lsst2016-i",
        }
        ### LSST expectations baked in
        kwargs_quasar = {
            "variability_model": "light_curve",
            "kwargs_variability": {"agn_lightcurve", "u", "g", "r", "i", "z", "y"},
            "agn_driving_variability_model": "bending_power_law_from_distribution",
            "agn_driving_kwargs_variability": variable_agn_kwarg_dict,
            "lightcurve_time": np.linspace(0, length_of_light_curve, length_of_light_curve),
            "corona_height": 10,
            "r_resolution": 500,
        }
        self.source_quasar_plus_galaxies = sources.PointPlusExtendedSources(
            point_plus_extended_sources_list=self.quasar_source_plus_galaxy,
            cosmo=self.cosmo,
            sky_area=self.quasar_sky_area,
            kwargs_cut=self.kwargs_source_cut,
            list_type="astropy_table",
            catalog_type="skypy",
            point_source_type="quasar",
            extended_source_type="single_sersic",
            point_source_kwargs=kwargs_quasar,
        )
    def generate_lens_population(self, image_separation=0.5, magnitude_limit=24, n_iter=1):
    # Initiate LensPop class to generate lensed quasar pop.
        quasar_lens_pop_ell = LensPop(
            deflector_population=self.lens_galaxies_ell,
            source_population=self.source_quasar_plus_galaxies,
            cosmo=self.cosmo,
            sky_area=self.full_sky_area,
        )
        kwargs_lens_cuts = {
            "min_image_separation": image_separation,
            "max_image_separation": 10,
            "second_brightest_image_cut": {"i": magnitude_limit},
        }
        # drawing population
        # the key difference in lens population drawing time is whether you ask for magnitude cuts or not I think?
        quasar_lens_population = []
        for _ in tqdm(range(n_iter)):
            qlp5000 = quasar_lens_pop_ell.draw_population(
                speed_factor=1000, kwargs_lens_cuts=kwargs_lens_cuts
            )
            quasar_lens_population.extend(qlp5000)

        return quasar_lens_population
    
    def generate_lens_population_dataset(self, quasar_lens_pop):
        f = io.StringIO()
        quasar_lens_pop_copy = copy.deepcopy(quasar_lens_pop)
        full_pop_df = pd.DataFrame()
        with redirect_stdout(f):
            for i, lens_obj in tqdm(enumerate(quasar_lens_pop)):
                full_pop_df = lens_obj.lens_to_dataframe(index=i, df=full_pop_df)
                image2mag = full_pop_df.loc[i, "point_source_light_i_magnitude_1"]
                try:
                    image3mag = full_pop_df.loc[i, "point_source_light_i_magnitude_2"]
                except KeyError:
                    image3mag = 0
                second_or_third_mag = (
                    image3mag if not (np.isnan(image3mag) or image3mag == 0) else image2mag)
                
                full_pop_df.loc[i, "i2"] = image2mag
                full_pop_df.loc[i, "i3"] = second_or_third_mag
                (
                    full_pop_df.loc[i, "deflector_mass_phi"],
                    full_pop_df.loc[i, "deflector_mass_q"],
                ) = ellipticity2phi_q(
                    full_pop_df.loc[i, "deflector_mass_e1"],
                    full_pop_df.loc[i, "deflector_mass_e2"],
                )
                (
                    full_pop_df.loc[i, "deflector_light_phi"],
                    full_pop_df.loc[i, "deflector_light_q"],
                ) = ellipticity2phi_q(
                    full_pop_df.loc[i, "deflector_light_i_e1"],
                    full_pop_df.loc[i, "deflector_light_i_e2"],
                )
                full_pop_df.loc[i, "deflector_stellar_mass"] = lens_obj.deflector_stellar_mass()
                full_pop_df.loc[i, "lens_obj"] = lens_obj
                for band in list("ugrizy"):
                    abs_mag = self.quasar_class.convert_magnitude(
                        full_pop_df.loc[i, f"ps_{band}_mag_true"],
                        full_pop_df.loc[i, "point_source_redshift"],
                        conversion="apparent_to_absolute",
                    )
                    full_pop_df.loc[i, f"M_{band}"] = abs_mag
                    if np.isnan(abs_mag):
                        if band == "y":
                            full_pop_df.drop(index=i, inplace=True)
                            quasar_lens_pop_copy.pop(i)
        max_num_of_images_in_df = int(np.max(full_pop_df['num_ps_images']))
        mask = (np.array(full_pop_df[[f"micro_kappa_star_{i}" for i in range(max_num_of_images_in_df)]]) 
                >= np.array(full_pop_df[[f"micro_kappa_tot_{i}" for i in range(max_num_of_images_in_df)]])).any(axis=1)
        full_pop_df = full_pop_df.loc[~mask].reset_index(drop = True)
        quasar_lens_pop_copy = np.array(quasar_lens_pop_copy)[~mask]
        self.lens_population_dataframe = full_pop_df
        self.quasar_lens_population = quasar_lens_pop_copy

    def generate_and_save_single_lens(self, lens_index, baseline=10, filename='lens_finding_postage_stamps.h5'):
        """
        Generate all necessary data for a lens and save it to HDF5.

        Parameters:
        -----------
        lens_index : int
            Index of the lens to generate
        baseline : int
            Baseline in years for observations (default: 10)
        filename : str
            Path to HDF5 file (default: 'lens_finding_postage_stamps.h5')
        """
        
        # Get metadata from dataframe
        lens_row = self.lens_population_dataframe.iloc[lens_index]
        
        # Get lens object and metadata from index
        lens_obj = self.quasar_lens_population[lens_index]    
        # Get random RA and Dec for observation
        new_ra, new_dec = get_random_ra_dec(N=1)

        
        bands = list("ugrizy")
        mag_zps = np.array(
            [
                26.52,
                28.51,
                28.36,
                28.17,
                27.78,
                26.82,
            ]
        )  # taken from https://smtn-002.lsst.io/

        mag_zero_points_30_seconds = dict(
            zip(bands, compute_magnitude_zeropoint(mag_zps))
        )  # mag
        delta_pix = 0.2  # arcsec/pixel
        num_pix = 33  # pixels
        exp_time = 30  # s
        
        # Get observation dates
        if baseline == 10:
            sql = ''
        else:
            sql = f'night < {baseline*365.25}'
        while True:
            try:    
                rubin_df = get_rubin_cadence(new_ra, new_dec, sql=sql)
                break
            except Exception as e:
                new_ra, new_dec = get_random_ra_dec(N=1)
                print(f"get_rubin_cadence failed, retrying: {e}")
        self.lens_population_dataframe.loc[lens_index, 'RA'] = new_ra.value
        self.lens_population_dataframe.loc[lens_index, 'DEC'] = new_dec.value
        observation_dates = rubin_df["observationStartMJD"]
        image_number = lens_obj.image_number
        
        if isinstance(image_number, list):
            image_number = image_number[0]
        # Generate postage stamp images for all bands
        image_lens_series_all_bands = []
        light_curves_dict = {}
        for img_idx in range(image_number):
            light_curves_dict[img_idx] = {}
        image_static_all_bands = []
        for band in bands:
            try:
                time_sampled = np.array(observation_dates[band])
                repeats = len(time_sampled)
                transform_matrix = np.array([[delta_pix, 0], [0, delta_pix]])
                psf_kernel_list = [None]
                transform_matrix_list = [transform_matrix]
                mag_list = [mag_zero_points_30_seconds[band]]
                expo_list = [exp_time]
                mag_zero_points_all = mag_list * repeats
                psf_kernels_all = psf_kernel_list * repeats
                transform_matrix_all = transform_matrix_list * repeats
                exposure_time_all = expo_list * repeats
                static_image = lens_image(
                    lens_class=lens_obj,
                    band=band,
                    mag_zero_point=mag_zero_points_30_seconds[band],
                    num_pix=num_pix,
                    psf_kernel=None,
                    transform_pix2angle=np.array([[delta_pix, 0], [0, delta_pix]]),
                    exposure_time=exp_time,
                    t_obs=None,
                    std_gaussian_noise=None,
                    with_source=True,
                    with_ps=False,
                    with_deflector=True,
                    add_noise=False,
                    gain=0.7,
                    single_visit_mag_zero_points={
                        "g": 32.33,
                        "r": 32.17,
                        "i": 31.85,
                        "z": 31.45,
                        "y": 30.63,
                    },
                )
                image_static_all_bands.append(static_image)
                image_lens_series, magnitudes = lens_image_series_precomputed_mags(
                    lens_class=lens_obj,
                    band=band,
                    mag_zero_point=mag_zero_points_all,
                    num_pix=num_pix,
                    psf_kernel=psf_kernels_all,
                    transform_pix2angle=transform_matrix_all,
                    exposure_time=exposure_time_all,
                    std_gaussian_noise=bkg_noise(
                        0.005, 30, np.array(rubin_df.loc[band, "skyBrightness"]), 0.2, 1
                    ),
                    t_obs=time_sampled,
                    with_deflector=True,
                    with_ps=True,
                    with_source=True,
                    add_noise=False,
                    single_visit_mag_zero_points=mag_zero_points_30_seconds
                )
                magnitudes = np.array(magnitudes).T
                image_lens_series_all_bands.append(image_lens_series)
                for img_idx in range(image_number):
                    light_curves_dict[img_idx][band] = magnitudes[img_idx]
            except KeyError:
                print(f'no observation {band} band in {baseline} years')
        metadata_dict = lens_row.to_dict()
        save_lens_to_hdf5(
            filename=filename,
            lens_index=lens_index, 
            metadata_dict=metadata_dict,
            observation_dates=observation_dates,
            static_images_all_bands=image_static_all_bands,
            image_lens_series_all_bands=image_lens_series_all_bands,
            light_curves_dict=light_curves_dict)
        
        return lens_obj
    
def save_lens_to_hdf5(
    filename,
    lens_index,
    metadata_dict,
    observation_dates,
    static_images_all_bands,
    image_lens_series_all_bands,
    light_curves_dict,
):
    """
    Save lens system data to HDF5 file.
    
    Parameters:
    -----------
    filename : str
        Path to HDF5 file (will be created if doesn't exist)
    lens_index : int
        Index of the lens (e.g., 1 for lsst_lens_1)
    metadata_dict : dict
        Dictionary containing metadata from full_pop_df row
    observation_dates : dict
        Dictionary with keys as bands ('u', 'g', 'r', 'i', 'z', 'y') and values as observation MJD arrays
    image_lens_series_all_bands : list
        List of image arrays for each band (6 bands total)
    light_curves_dict : dict
        Dictionary with structure {image_num: {band: magnitude_array}}
        e.g., {0: {'u': [mag1, mag2, ...], 'g': [mag1, mag2, ...]}, 1: {...}, ...}
    """
    
    # Open file in append mode (creates if doesn't exist)
    with h5py.File(filename, 'a') as hf:
        # Create group for this lens
        lens_group_name = f'lsst_lens_{lens_index}'
        
        # Delete group if it already exists
        if lens_group_name in hf:
            del hf[lens_group_name]
        
        lens_group = hf.create_group(lens_group_name)
        
        # Save metadata
        metadata_group = lens_group.create_group('metadata')
        for key, value in metadata_dict.items():
            # Handle different data types
            if isinstance(value, (int, float, np.integer, np.floating)):
                metadata_group.attrs[key] = value
            elif isinstance(value, str):
                metadata_group.attrs[key] = value
            elif isinstance(value, np.ndarray):
                metadata_group.create_dataset(key, data=value)
            elif value is None or (isinstance(value, float) and np.isnan(value)):
                metadata_group.attrs[key] = 'NaN'
            else:
                # Try to convert to string as fallback
                try:
                    metadata_group.attrs[key] = str(value)
                except:
                    print(f"Warning: Could not save metadata key '{key}' with value {value}")
        
        # Save observation dates for each band
        obs_dates_group = lens_group.create_group('observation_dates')
        bands = ['u', 'g', 'r', 'i', 'z', 'y']
        for band in bands:
            if band in observation_dates:
                obs_dates_group.create_dataset(band, data=np.array(observation_dates[band]))
        
        # Save postage stamp images for each band
        images_group = lens_group.create_group('postage_stamps')
        for i, band in enumerate(bands):
            if i < len(image_lens_series_all_bands):
                band_group = images_group.create_group(band)
                image_series = np.array(image_lens_series_all_bands[i])
                band_group.create_dataset(f'all_observations', data=image_series, compression='gzip')
        static_images_group = lens_group.create_group('static_image')
        for i, band in enumerate(bands):
            if i < len(static_images_all_bands):
                band_group = static_images_group.create_group(band)
                image_band = np.array(static_images_all_bands[i])
                band_group.create_dataset(f'lens_plus_lensed_agn_host', data=image_band, compression='gzip')
        
        # Save light curves for each band and image
        lightcurves_group = lens_group.create_group('light_curves')
        # Save light curves for each image and band
        for image_num, bands_lc in light_curves_dict.items():
            image_group = lightcurves_group.create_group(f'image_{image_num}')
            for band, magnitudes in bands_lc.items():
                image_group.create_dataset(band, data=np.array(magnitudes))

def main():
    """
    Main function to generate and save lensed quasar-host galaxy systems for LSST.
    Handles command-line arguments and orchestrates the full dataset generation pipeline.
    """
    parser = argparse.ArgumentParser(
        description='Generate and save LSST lensed quasar postage stamp dataset'
    )
    parser.add_argument(
        '--sky-area',
        type=float,
        default=10,
        help='Sky area in deg^2 for lens sampling (default: 10)'
    )
    parser.add_argument(
        '--n-iter',
        type=int,
        default=1,
        help='Number of iterations for lens population generation (default: 1)'
    )
    parser.add_argument(
        '--image-separation',
        type=float,
        default=0.5,
        help='Minimum image separation in arcsec (default: 0.5)'
    )
    parser.add_argument(
        '--magnitude-limit',
        type=float,
        default=24,
        help='Magnitude limit for second brightest image in i-band (default: 24)'
    )
    parser.add_argument(
        '--baseline',
        type=int,
        default=10,
        help='Baseline in years for observations (default: 10)'
    )
    parser.add_argument(
        '--output-file',
        type=str,
        default='lens_finding_postage_stamps.h5',
        help='Output HDF5 filename (default: lens_finding_postage_stamps.h5)'
    )
    parser.add_argument(
        '--start-index',
        type=int,
        default=0,
        help='Starting lens index (default: 0)'
    )
    parser.add_argument(
        '--end-index',
        type=int,
        default=None,
        help='Ending lens index (optional, if not provided, uses all generated lenses)'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        default='lens_generation_log.json',
        help='Log file to track failed lenses and generation progress (default: lens_generation_log.json)'
    )
    
    args = parser.parse_args()
    
    print(f"Starting LSST lensed quasar dataset generation")
    print(f"Configuration:")
    print(f"  Sky area: {args.sky_area} deg^2")
    print(f"  N iterations: {args.n_iter}")
    print(f"  Image separation: {args.image_separation} arcsec")
    print(f"  Magnitude limit: {args.magnitude_limit}")
    print(f"  Baseline: {args.baseline} years")
    print(f"  Output file: {args.output_file}")
    print()
    
    # Initialize dataset generator
    generator = DatasetGenerator(
        sky_area=args.sky_area,
        bands=list('ugrizy'),
        postage_stamps_path='.',
        dataset_path='.'
    )
    
    # Set up cosmology
    print("Setting up cosmology...")
    generator.set_up_cosmology()
    
    # Set up galaxy and quasar populations
    print("Generating galaxy and quasar populations...")
    generator.set_up_galaxy_quasar_population()
    
    # Generate lens population
    print(f"Generating lens population with {args.n_iter} iteration(s)...")
    quasar_lens_pop = generator.generate_lens_population(
        image_separation=args.image_separation,
        magnitude_limit=args.magnitude_limit,
        n_iter=args.n_iter
    )
    print(f"Generated {len(quasar_lens_pop)} lens systems")
    
    # Generate lens population dataset (metadata)
    print("Computing lens population metadata...")
    generator.generate_lens_population_dataset(quasar_lens_pop)
    metadata_csv = args.output_file.split('.')[0] + '.csv'
    generator.lens_population_dataframe.to_csv(metadata_csv)
    # Determine lens indices to process
    num_lenses = len(generator.quasar_lens_population)
    end_index = args.end_index if args.end_index is not None else num_lenses
    end_index = min(end_index, num_lenses)
    
    lens_indices = range(args.start_index, end_index)
    print(f"Processing lenses {args.start_index} to {end_index-1} (total: {len(lens_indices)})")
    print()
    
    # Initialize logging
    log_data = {
        'timestamp': datetime.now().isoformat(),
        'config': vars(args),
        'successful_lenses': [],
        'failed_lenses': [],
        'total_processed': 0,
    }
    
    # Generate and save each lens with error handling
    for lens_idx in tqdm(lens_indices, desc="Processing lenses"):
        try:
            generator.generate_and_save_single_lens(
                lens_index=lens_idx,
                baseline=args.baseline,
                filename=args.output_file
            )
            log_data['successful_lenses'].append(lens_idx)
            print(f"  ✓ Lens {lens_idx} saved successfully")
        except Exception as e:
            log_data['failed_lenses'].append({
                'lens_index': lens_idx,
                'error': str(e),
                'error_type': type(e).__name__
            })
            print(f"  ✗ Lens {lens_idx} failed: {type(e).__name__}: {e}")
        finally:
            log_data['total_processed'] += 1
    
    # Save log file
    print()
    print(f"Saving log to {args.log_file}...")
    with open(args.log_file, 'w') as f:
        json.dump(log_data, f, indent=2)
    
    # Print summary
    print()
    print("=" * 60)
    print("GENERATION SUMMARY")
    print("=" * 60)
    print(f"Total processed: {log_data['total_processed']}")
    print(f"Successful: {len(log_data['successful_lenses'])}")
    print(f"Failed: {len(log_data['failed_lenses'])}")
    if log_data['failed_lenses']:
        print()
        print("Failed lens indices:")
        for failure_info in log_data['failed_lenses']:
            print(f"  - Lens {failure_info['lens_index']}: {failure_info['error_type']}")
    print("=" * 60)


if __name__ == "__main__":
    main()
