#!/usr/bin/env python
# coding: utf-8

# 
# # Lantern individual-target HTML viewer — v21 test
# 
# This notebook (previous name `lantern_individual_target_html_viewer_v21`) generates **one self-contained HTML page for one Second-Pass `target_id`**.
# 
# Main candidate input:
# 
# ```text
# candidate_analysis_hybrid_multi_loci_v8_no_gaia_color_cut/
#     second_pass_candidates_with_flags.parquet
# ```
# 
# Main features:
# 
# - target metadata from the Gaia-cleaning candidate table;
# - Galactic latitude;
# - First-Cut and Second-Pass scores;
# - Gaia PM S/N, parallax S/N, BP-RP, GaiaQSO, Milliquas, KnownLens;
# - `ugrizy` magnitudes from the **most recent alert's**
#   `lsst_diaObject_{band}_scienceFluxMean`;
# - alert count after 2026-05-27 and per-band counts;
# - **median DiaSource extendedness computed directly from the target's alerts**;
# - DP2 IRG and ZIR deep-coadd color images;
# - Legacy Survey cutout;
# - per-band DP2 science cutouts + local PSF kernels for future deconvolution;
# - placeholder multiband deconvolution and deconvolved IRG/ZIR;
# - ANTARES **science/visit + template + difference** thumbnails for as many
#   post-2026-05-27 alerts as can be retrieved;
# - a synchronized JavaScript animation of the three thumbnail types.
# 
# ### Important change from v4
# 
# The DP2 `Visit` table is **not used**. There is therefore:
# 
# - no `skyRotation` query;
# - no image de-rotation;
# - no TAP dependency for the viewer.
# 
# The visit number shown in the HTML comes only from the alert property
# `lsst_diaSource_visit`.
# 

# 
# ## Notes for v21
# 
# ### DiaSource extendedness
# 
# This notebook computes the displayed median directly from:
# 
# ```text
# lsst_diaSource_extendedness
# ```
# 
# in the saved ANTARES alert properties, using only the target's retained alerts
# on/after 2026-05-27.
# 
# The Rubin APDB/DP2 schema defines `extendedness` as a **DiaSource** quantity,
# not a DiaObject quantity. It compares the source's moment-based size to the PSF
# moments; values near 1 indicate an extended source and values near 0 indicate a
# point-like source.
# 
# ### ANTARES thumbnails
# 
# For every retained alert, the notebook calls:
# 
# ```python
# from antares_client.search import get_thumbnails
# get_thumbnails("lsst:<alert_id>")
# ```
# 
# and attempts to retrieve all three alert cutouts:
# 
# - science / visit
# - template
# - difference
# 
# Alerts are tried newest first. Failures are recorded and do not stop retrieval
# for later alerts.
# 
# ### Pixel display
# 
# Small LSST cutouts are shown with CSS nearest-neighbor / `pixelated` rendering
# so the original detector/coadd pixels remain visibly sharp when the browser
# enlarges them.
# 
# 
# ### Missing-band handling
# 
# The notebook now distinguishes between:
# 
# 1. bands listed by the DP2 SIA search; and
# 2. bands for which the actual exposure cutout **and** local PSF were successfully retrieved.
# 
# The final `requested -> source` mapping is built only from the successfully
# retrieved bands. Therefore, if (for example) `u` is listed by SIA but its
# cutout fails, the requested `u` plane falls back to the nearest usable band
# instead of stopping the notebook.
# 
# 
# ### HTML fallback display
# 
# The HTML explicitly reports when a requested band was replaced by another usable
# band. This appears both in a top summary metric and in the per-band
# deconvolution rows.
# 

# In[1]:


# ============================================================
# USER INPUT
# ============================================================

# Set the target ID here.
import sys
TARGET_ID = sys.argv[1] #"target_0002330"


# ============================================================
# MAIN INPUTS
# ============================================================

from pathlib import Path

CANDIDATES_FILE = Path(
    "candidate_analysis_hybrid_multi_loci_v8_no_gaia_color_cut/"
    "second_pass_candidates_with_flags.parquet"
)

CANDIDATE_ALERTS_FILE = Path(
    "second_pass_output_hybrid_multi_loci_v3/"
    "second_pass_candidate_alerts.parquet"
)




# Optional known-lens catalog.
KNOWN_LENS_FILE = Path(
    "lenses_final.csv"
)

KNOWN_RA_COL = "ra"
KNOWN_DEC_COL = "dec"
KNOWN_NAME_COL = None

KNOWN_LENS_MATCH_RADIUS_ARCSEC = 5.0


# ============================================================
# OUTPUT
# ============================================================

OUTPUT_DIR = Path(
    "individual_target_html_viewer"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# ALERT / VIEWER SETTINGS
# ============================================================

MJD_MIN = 61187.0       # 2026-05-27
# Try all retained target alerts by default.
# Set to an integer for a quick test, e.g. 10.
MAX_THUMBNAIL_ALERTS = None

# Print one retrieval-status line per alert.
THUMBNAIL_VERBOSE = True

# DP2 coadd radius:
# 0.001 deg = 3.6 arcsec radius = 7.2 arcsec full width.
COADD_CUTOUT_RADIUS_DEG = 0.001


# ============================================================
# DECONVOLUTION PREPARATION
# ============================================================

DECONV_REQUESTED_BANDS = list("ugrizy")

# Full DP2 Exposure cutouts are used here so the PSF model is included.
DECONV_CUTOUT_RADIUS_DEG = COADD_CUTOUT_RADIUS_DEG

# Save science image + normalized PSF kernel for future deconvolution.
SAVE_DECONV_INPUT_NPZ = True

# Display-only size for the PSF panel in the HTML.
PSF_DISPLAY_PIX = 180

# ============================================================
# DECONVOLUTION PLACEHOLDER CONTROLS
# ============================================================

# Turn this on if you want the notebook to run the placeholder
# multiband deconvolution function and build deconvolved IRG/ZIR.
RUN_DECONVOLUTION = True

# Increase/decrease fake-color saturation for IRG/ZIR display.
# 1.0 = unchanged; 1.3–1.8 gives stronger color.
RGB_SATURATION = 1.6


# Legacy Survey image.
LEGACY_LAYER = "ls-dr11"
LEGACY_SIZE_PIX = 64
LEGACY_PIXSCALE = 0.262


SHOW_ALERT_ID = True

# JavaScript frame time.
ANIMATION_INTERVAL_MS = 850


# ============================================================
# ANTARES API
# ============================================================

ANTARES_BASE_URL = "https://api.antares.noirlab.edu/v2/"
ANTARES_ALERT_LIMIT = 1000

REQUEST_TIMEOUT = 60


# In[2]:


# ============================================================
# IMPORTS
# ============================================================

import ast
import base64
import html
import io
import json
import re
from collections import Counter
from collections.abc import Mapping
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests

from PIL import Image, ImageOps, ImageEnhance

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import AsinhStretch, LinearStretch, ImageNormalize

from IPython.display import display, FileLink, Markdown
from antares_client.search import get_thumbnails

from lsst.images.serialization import read_archive
from lsst.rsp import RSPDiscovery
from lsst.rsp.utils import get_pyvo_auth
from pyvo.dal.adhoc import SodaQuery


# In[3]:


sys.path.append('/home/sfu/SL_time_delay_2026a/loci/SharPy')
from sharpy.enhancement import Enhancement

my_enh = Enhancement(device='GPU', trained_model_file=None)


# ## 1. Load the Gaia-cleaning candidate table and select one target

# In[4]:


candidates = pd.read_parquet(
    CANDIDATES_FILE
)

if "target_id" not in candidates.columns:
    raise KeyError(
        "The Gaia-cleaning candidate table must contain target_id."
    )

matches = candidates[
    candidates["target_id"].astype(str)
    == str(TARGET_ID)
].copy()

if len(matches) == 0:
    raise KeyError(
        f"TARGET_ID={TARGET_ID!r} was not found in {CANDIDATES_FILE}."
    )

if len(matches) > 1:
    raise RuntimeError(
        f"TARGET_ID={TARGET_ID!r} appears {len(matches)} times."
    )

target = matches.iloc[0].copy()

print(
    f"Selected target: {target['target_id']}"
)

print(
    f"RA, Dec: {target['ra']:.8f}, {target['dec']:.8f}"
)


# ## 2. Load the target's alerts after 2026-05-27

# In[5]:


candidate_alerts = pd.read_parquet(
    CANDIDATE_ALERTS_FILE
)

required_alert_columns = {
    "target_id",
    "locus_id",
    "alert_id",
    "mjd",
    "alert_properties",
}

missing = (
    required_alert_columns
    - set(candidate_alerts.columns)
)

if missing:
    raise KeyError(
        "Candidate-alert file is missing required columns: "
        f"{sorted(missing)}"
    )


def parse_properties(value):
    if isinstance(value, Mapping):
        return dict(value)

    if isinstance(value, str):
        for parser in (
            json.loads,
            ast.literal_eval,
        ):
            try:
                parsed = parser(value)
                if isinstance(parsed, Mapping):
                    return dict(parsed)
            except Exception:
                pass

    return {}


target_alerts = candidate_alerts[
    candidate_alerts["target_id"].astype(str)
    == str(TARGET_ID)
].copy()

target_alerts["mjd"] = pd.to_numeric(
    target_alerts["mjd"],
    errors="coerce",
)

target_alerts = target_alerts[
    target_alerts["mjd"] >= MJD_MIN
].copy()

target_alerts["_props"] = target_alerts[
    "alert_properties"
].map(
    parse_properties
)


def get_prop(props, keys, default=np.nan):
    for key in keys:
        if key in props:
            return props[key]
    return default


def to_float(value):
    try:
        return float(value)
    except Exception:
        return np.nan


def to_int_or_none(value):
    try:
        if value is None:
            return None

        x = int(value)
        return x

    except Exception:
        try:
            x = int(float(value))
            return x
        except Exception:
            return None


target_alerts["band"] = target_alerts[
    "_props"
].map(
    lambda p: get_prop(
        p,
        [
            "lsst_diaSource_band",
            "diaSource_band",
            "band",
        ],
        default=None,
    )
)

target_alerts["visit"] = target_alerts[
    "_props"
].map(
    lambda p: to_int_or_none(
        get_prop(
            p,
            [
                "lsst_diaSource_visit",
                "diaSource_visit",
                "visit",
            ],
            default=None,
        )
    )
)

target_alerts = target_alerts.sort_values(
    "mjd",
    ascending=False,
).reset_index(
    drop=True
)

if len(target_alerts) == 0:
    raise RuntimeError(
        "No target alerts remain after the 2026-05-27 MJD cut."
    )

latest_alert = target_alerts.iloc[0]

print(
    f"N alerts after 2026-05-27: {len(target_alerts)}"
)

print(
    "Latest alert:"
)

display(
    target_alerts[
        [
            "locus_id",
            "alert_id",
            "mjd",
            "band",
            "visit",
        ]
    ].head(
        10
    )
)


# 
# ## 3. Compute target summary quantities
# 
# `extendedness_median` is computed directly from the retained alert property
# 
# ```text
# lsst_diaSource_extendedness
# ```
# 
# for this target. No separate diagnostics output directory is required.
# 

# In[6]:


ra = float(
    target["ra"]
)

dec = float(
    target["dec"]
)

coord = SkyCoord(
    ra=ra * u.deg,
    dec=dec * u.deg,
    frame="icrs",
)

galactic_latitude_deg = float(
    coord.galactic.b.deg
)

latest_locus_id = str(
    latest_alert["locus_id"]
)


# ------------------------------------------------------------
# DiaSource extendedness from retained alerts.
# ------------------------------------------------------------

target_alerts["diaSource_extendedness"] = target_alerts[
    "_props"
].map(
    lambda p: to_float(
        p.get(
            "lsst_diaSource_extendedness",
            np.nan,
        )
    )
)

finite_extendedness = target_alerts.loc[
    np.isfinite(
        target_alerts[
            "diaSource_extendedness"
        ]
    ),
    "diaSource_extendedness",
].to_numpy(
    dtype=float
)

extendedness_n = len(
    finite_extendedness
)

extendedness_median = (
    float(
        np.median(
            finite_extendedness
        )
    )
    if extendedness_n > 0
    else np.nan
)


# ------------------------------------------------------------
# ugrizy magnitudes from the MOST RECENT alert only.
# ------------------------------------------------------------

AB_ZEROPOINT_NJY = 31.4


def flux_njy_to_abmag(flux):
    flux = to_float(
        flux
    )

    if (
        not np.isfinite(flux)
        or flux <= 0
    ):
        return np.nan

    return (
        AB_ZEROPOINT_NJY
        - 2.5
        * np.log10(
            flux
        )
    )


latest_props = latest_alert[
    "_props"
]

latest_alert_mags = {}

for band in "ugrizy":

    flux = to_float(
        latest_props.get(
            f"lsst_diaObject_{band}_scienceFluxMean",
            np.nan,
        )
    )

    latest_alert_mags[
        band
    ] = flux_njy_to_abmag(
        flux
    )

print(
    ",".join(
        [
            str(TARGET_ID).split('_')[1],
            f"{ra:.8f}",
            f"{dec:.8f}",
        ]
        + [
            (
                f"{latest_alert_mags[band]:.4f}"
                if np.isfinite(latest_alert_mags[band])
                #else "nan"
                else "--"
            )
            for band in "ugrizy"
        ]
    )
)

# ------------------------------------------------------------
# Band counts after May 27.
# ------------------------------------------------------------

band_counts = Counter(
    str(x)
    for x
    in target_alerts["band"]
    if x is not None
    and str(x) != "nan"
)


print(
    f"Galactic latitude b = {galactic_latitude_deg:.4f} deg"
)

print(
    f"Latest-alert locus ID = {latest_locus_id}"
)

print(
    "Band counts:",
    dict(
        band_counts
    ),
)

print(
    "Most-recent-alert magnitudes:",
    latest_alert_mags,
)


print(
    "DiaSource extendedness:",
    f"N={extendedness_n}, median={extendedness_median}"
)


# ## 4. Optional KnownLens cross-match

# In[7]:


known_lens_match = None
known_lens_name = None
known_lens_sep_arcsec = np.nan


if KNOWN_LENS_FILE.exists():

    known = pd.read_csv(
        KNOWN_LENS_FILE
    )

    if (
        KNOWN_RA_COL in known.columns
        and KNOWN_DEC_COL in known.columns
    ):

        known = known.copy()

        known[KNOWN_RA_COL] = pd.to_numeric(
            known[KNOWN_RA_COL],
            errors="coerce",
        )

        known[KNOWN_DEC_COL] = pd.to_numeric(
            known[KNOWN_DEC_COL],
            errors="coerce",
        )

        known = known.dropna(
            subset=[
                KNOWN_RA_COL,
                KNOWN_DEC_COL,
            ]
        ).reset_index(
            drop=True
        )


        if len(
            known
        ):

            known_coord = SkyCoord(
                ra=known[KNOWN_RA_COL].to_numpy() * u.deg,
                dec=known[KNOWN_DEC_COL].to_numpy() * u.deg,
            )

            idx, sep2d, _ = coord.match_to_catalog_sky(
                known_coord
            )

            known_lens_sep_arcsec = float(
                sep2d.arcsec
            )

            known_lens_match = (
                known_lens_sep_arcsec
                <= KNOWN_LENS_MATCH_RADIUS_ARCSEC
            )

            if (
                known_lens_match
                and KNOWN_NAME_COL is not None
                and KNOWN_NAME_COL in known.columns
            ):
                known_lens_name = str(
                    known.iloc[int(idx)][KNOWN_NAME_COL]
                )


if known_lens_match is None:
    print(
        "KnownLens: N/A "
        "(known-lens catalog unavailable or missing RA/Dec columns)"
    )
else:
    print(
        f"KnownLens: {known_lens_match}, "
        f"nearest separation={known_lens_sep_arcsec:.3f} arcsec"
    )


# ## 5. Initialize DP2 SIA client

# In[8]:


discovery = RSPDiscovery(
    "dp2"
)

sia_client = discovery.get_sia_client()


# 
# ## 6. DP2 deep-coadd exposures, IRG/ZIR, and per-band PSF models
# 
# This section makes one SIA search at the target position, determines which LSST
# bands are available, and retrieves full exposure cutouts with
# `cutout-sync-exposure`.
# 
# For every requested band in `ugrizy`:
# 
# - use the exact band when available;
# - otherwise choose the nearest available band in filter order;
# - for equal-distance choices, prefer the band closer to the middle of `ugrizy`.
# 
# Examples:
# 
# ```text
# u missing -> g
# y missing -> z
# ```
# 
# The HTML explicitly labels every substitution.
# 

# In[9]:


BAND_ORDER = list("ugrizy")
BAND_INDEX = {
    band: i
    for i, band in enumerate(BAND_ORDER)
}
BAND_CENTER_INDEX = (len(BAND_ORDER) - 1) / 2.0


# ------------------------------------------------------------
# One SIA search for all deep-coadd bands at the target.
# ------------------------------------------------------------

coadd_search = sia_client.search(
    pos=(
        ra,
        dec,
        COADD_CUTOUT_RADIUS_DEG,
    ),
    calib_level=3,
    dpsubtype="lsst.deep_coadd",
)

coadd_table = coadd_search.to_table()

available_bands = sorted(
    {
        str(x)
        for x in coadd_table["lsst_band"]
        if str(x) in BAND_ORDER
    },
    key=lambda b: BAND_INDEX[b],
)

if not available_bands:
    raise RuntimeError(
        "No DP2 deep_coadd bands were found at this target."
    )

print(
    "Available DP2 deep-coadd bands:",
    available_bands,
)


def choose_source_band(
    requested_band,
    available_bands,
):
    """
    Use exact band if available.

    Otherwise choose the nearest band in ugrizy order.
    For an equal-distance tie, prefer the band nearer the
    middle of ugrizy.
    """
    requested_band = str(
        requested_band
    )

    if requested_band in available_bands:
        return requested_band

    req_i = BAND_INDEX[
        requested_band
    ]

    ranked = sorted(
        available_bands,
        key=lambda band: (
            abs(
                BAND_INDEX[band]
                - req_i
            ),
            abs(
                BAND_INDEX[band]
                - BAND_CENTER_INDEX
            ),
            BAND_INDEX[band],
        ),
    )

    return ranked[0]



# Final requested->source mapping is determined only AFTER
# the actual cutout + PSF retrieval succeeds.

def best_coadd_result_index_for_band(
    band,
):
    indices = np.where(
        np.asarray(
            coadd_table["lsst_band"]
        ).astype(str)
        == str(band)
    )[0]

    if len(indices) == 0:
        return None

    distances = []

    for idx in indices:
        try:
            s_ra = float(
                coadd_table["s_ra"][idx]
            )

            s_dec = float(
                coadd_table["s_dec"][idx]
            )

            distance = (
                (s_ra - ra)**2
                + (s_dec - dec)**2
            )**0.5

        except Exception:
            distance = np.inf

        distances.append(
            distance
        )

    return int(
        indices[
            int(
                np.argmin(
                    distances
                )
            )
        ]
    )


EXPOSURE_CACHE = {}


def get_dp2_exposure_cutout(
    band,
    radius=DECONV_CUTOUT_RADIUS_DEG,
):
    band = str(
        band
    )

    if band in EXPOSURE_CACHE:
        return EXPOSURE_CACHE[
            band
        ]

    indices = np.where(
        np.asarray(
            coadd_table["lsst_band"]
        ).astype(str)
        == band
    )[0]

    if len(indices) == 0:
        EXPOSURE_CACHE[
            band
        ] = None
        return None

    best_idx = best_coadd_result_index_for_band(
        band
    )

    candidate_indices = []

    if best_idx is not None:
        candidate_indices.append(
            int(best_idx)
        )

    candidate_indices.extend(
        int(idx)
        for idx in indices
        if int(idx) not in candidate_indices
    )

    errors = []

    for idx in candidate_indices:

        try:
            dl_result = discovery.get_datalink_results(
                coadd_search[
                    idx
                ]
            )

            sq = SodaQuery.from_resource(
                dl_result,
                dl_result.get_adhocservice_by_id(
                    "cutout-sync-exposure"
                ),
                session=get_pyvo_auth(),
            )

            sq.circle = (
                ra * u.deg,
                dec * u.deg,
                radius * u.deg,
            )

            cutout_bytes = (
                sq.execute_stream()
                .read()
            )

            sq.raise_if_error()

            exposure = read_archive(
                io.BytesIO(
                    cutout_bytes
                )
            )

            EXPOSURE_CACHE[
                band
            ] = exposure

            return exposure

        except Exception as exc:
            errors.append(
                f"coadd index {idx}: {exc}"
            )

    EXPOSURE_CACHE[
        band
    ] = None

    raise RuntimeError(
        f"All DP2 deep-coadd cutout attempts failed for band {band}. "
        + " | ".join(
            errors
        )
    )


def get_exposure_image_array(
    exposure
):
    if exposure is None:
        return None

    if hasattr(
        exposure,
        "image",
    ):
        return np.asarray(
            exposure.image.array,
            dtype=float,
        )

    if hasattr(
        exposure,
        "array",
    ):
        return np.asarray(
            exposure.array,
            dtype=float,
        )

    raise AttributeError(
        "Could not find image pixel array on DP2 exposure."
    )


def target_pixel_in_exposure(
    exposure
):
    try:
        target_coord = SkyCoord(
            ra=ra * u.deg,
            dec=dec * u.deg,
            frame="icrs",
        )

        xy = (
            exposure
            .sky_projection
            .sky_to_pixel(
                target_coord
            )
        )

        return (
            float(xy.x),
            float(xy.y),
        )

    except Exception:
        bbox = exposure.bbox

        x = (
            bbox.x.min
            + (
                bbox.x.max
                - bbox.x.min
            )
            / 2.0
        )

        y = (
            bbox.y.min
            + (
                bbox.y.max
                - bbox.y.min
            )
            / 2.0
        )

        return (
            float(x),
            float(y),
        )


def get_local_psf_kernel_array(
    exposure
):
    if exposure is None:
        return None

    x_pix, y_pix = target_pixel_in_exposure(
        exposure
    )

    psf_image = (
        exposure.psf
        .compute_kernel_image(
            x=x_pix,
            y=y_pix,
        )
    )

    psf_array = np.asarray(
        psf_image.array,
        dtype=float,
    )

    total = np.nansum(
        psf_array
    )

    if (
        np.isfinite(total)
        and total != 0
    ):
        psf_array = (
            psf_array
            / total
        )

    return psf_array



# ------------------------------------------------------------
# Try every band reported by SIA.
#
# Important: a band being present in the SIA table does not guarantee
# that the actual exposure cutout and PSF retrieval will succeed.
# ------------------------------------------------------------

band_products = {}
band_fetch_errors = {}


for source_band in available_bands:

    print(
        f"Fetching full DP2 exposure + PSF for {source_band}..."
    )

    try:

        exposure = get_dp2_exposure_cutout(
            source_band
        )

        image_array = get_exposure_image_array(
            exposure
        )

        psf_array = get_local_psf_kernel_array(
            exposure
        )

        if image_array is None:
            raise RuntimeError(
                "science image is None"
            )

        if psf_array is None:
            raise RuntimeError(
                "PSF array is None"
            )

        image_array = np.asarray(
            image_array,
            dtype=float,
        )

        psf_array = np.asarray(
            psf_array,
            dtype=float,
        )

        if image_array.ndim != 2:
            raise RuntimeError(
                f"science image is not 2D: shape={image_array.shape}"
            )

        if psf_array.ndim != 2:
            raise RuntimeError(
                f"PSF is not 2D: shape={psf_array.shape}"
            )

        if not np.any(
            np.isfinite(
                image_array
            )
        ):
            raise RuntimeError(
                "science image has no finite pixels"
            )

        if not np.any(
            np.isfinite(
                psf_array
            )
        ):
            raise RuntimeError(
                "PSF has no finite pixels"
            )

        band_products[
            source_band
        ] = {
            "exposure":
                exposure,
            "image":
                image_array,
            "psf":
                psf_array,
        }

        print(
            "  success:",
            "image shape",
            image_array.shape,
            "| PSF shape",
            psf_array.shape,
        )

    except Exception as exc:

        band_fetch_errors[
            source_band
        ] = str(
            exc
        )

        print(
            f"  failed for {source_band}: {exc}"
        )

        band_products[
            source_band
        ] = {
            "exposure":
                None,
            "image":
                None,
            "psf":
                None,
        }


# ------------------------------------------------------------
# Keep only bands for which BOTH science image and PSF worked.
# ------------------------------------------------------------

usable_bands = [
    band
    for band
    in BAND_ORDER
    if (
        band in band_products
        and band_products[
            band
        ][
            "image"
        ]
        is not None
        and band_products[
            band
        ][
            "psf"
        ]
        is not None
    )
]


print()
print(
    "Usable bands with both science image and PSF:",
    usable_bands,
)


if band_fetch_errors:

    print()
    print(
        "Bands that failed actual cutout/PSF retrieval:"
    )

    for band in BAND_ORDER:

        if band in band_fetch_errors:

            print(
                f"  {band}: {band_fetch_errors[band]}"
            )


if not usable_bands:

    raise RuntimeError(
        "No DP2 band has both a usable science image and PSF. "
        "Cannot construct the six-band deconvolution input."
    )


# ------------------------------------------------------------
# NOW build the requested -> source mapping from usable bands.
#
# Exact band is used when possible. Missing/failed bands fall back
# to the nearest usable band in ugrizy order.
# ------------------------------------------------------------

requested_to_source_band = {
    requested:
        choose_source_band(
            requested,
            usable_bands,
        )

    for requested
    in DECONV_REQUESTED_BANDS
}


print()
print(
    "Final requested -> usable source-band mapping:"
)

for requested in DECONV_REQUESTED_BANDS:

    source = requested_to_source_band[
        requested
    ]

    suffix = (
        ""
        if source == requested
        else "  [fallback]"
    )

    print(
        f"  {requested} -> {source}{suffix}"
    )

def grayscale_array_to_png_bytes(
    array,
):
    if array is None:
        return None

    data = np.asarray(
        array,
        dtype=float,
    )

    finite = data[
        np.isfinite(
            data
        )
    ]

    if len(finite) == 0:
        return None

    norm = ImageNormalize(
        data,
        stretch=LinearStretch(),
        clip=True,
    )

    scaled = np.asarray(
        norm(data),
        dtype=float,
    )

    scaled = np.nan_to_num(
        scaled,
        nan=0.0,
        posinf=1.0,
        neginf=0.0,
    )

    # Display-only orientation fix.
    scaled = np.flipud(
        scaled
    )

    image = Image.fromarray(
        np.uint8(
            np.round(
                255
                * np.clip(
                    scaled,
                    0,
                    1,
                )
            )
        ),
        mode="L",
    ).convert(
        "RGB"
    )

    buffer = io.BytesIO()

    image.save(
        buffer,
        format="PNG",
    )

    return buffer.getvalue()

def psf_array_to_png_bytes(
    psf_array,
):
    if psf_array is None:
        return None

    data = np.asarray(
        psf_array,
        dtype=float,
    )

    finite = data[
        np.isfinite(
            data
        )
    ]

    if len(finite) == 0:
        return None

    norm = ImageNormalize(
        data,
        stretch=LinearStretch(),
        clip=True,
    )

    scaled = np.asarray(
        norm(data),
        dtype=float,
    )

    scaled = np.nan_to_num(
        scaled,
        nan=0.0,
        posinf=1.0,
        neginf=0.0,
    )

    # Display-only orientation fix.
    scaled = np.flipud(
        scaled
    )

    image = Image.fromarray(
        np.uint8(
            np.round(
                255
                * np.clip(
                    scaled,
                    0,
                    1,
                )
            )
        ),
        mode="L",
    ).convert(
        "RGB"
    )

    buffer = io.BytesIO()

    image.save(
        buffer,
        format="PNG",
    )

    return buffer.getvalue()

deconv_band_rows = []


for requested_band in DECONV_REQUESTED_BANDS:

    source_band = requested_to_source_band[
        requested_band
    ]

    product = band_products.get(
        source_band,
        {},
    )

    image_array = product.get(
        "image"
    )

    psf_array = product.get(
        "psf"
    )

    deconv_band_rows.append(
        {
            "requested_band":
                requested_band,
            "source_band":
                source_band,
            "is_fallback":
                requested_band
                != source_band,
            "image":
                image_array,
            "psf":
                psf_array,
            "image_png":
                grayscale_array_to_png_bytes(
                    image_array
                ),
            "psf_png":
                psf_array_to_png_bytes(
                    psf_array
                ),
        }
    )


# ------------------------------------------------------------
# Save raw science image + normalized PSF arrays.
# ------------------------------------------------------------

target_id_safe_for_files = re.sub(
    r"[^A-Za-z0-9_.-]+",
    "_",
    str(TARGET_ID),
)

DECONV_INPUT_DIR = (
    OUTPUT_DIR
    / (
        f"target_{target_id_safe_for_files}"
        "_deconv_inputs"
    )
)

if SAVE_DECONV_INPUT_NPZ:

    DECONV_INPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    for row in deconv_band_rows:

        if (
            row["image"] is None
            or row["psf"] is None
        ):
            continue

        npz_path = (
            DECONV_INPUT_DIR
            / (
                f"{row['requested_band']}_"
                f"using_{row['source_band']}.npz"
            )
        )

        np.savez_compressed(
            npz_path,
            image=np.asarray(
                row["image"],
                dtype=float,
            ),
            psf=np.asarray(
                row["psf"],
                dtype=float,
            ),
            requested_band=np.asarray(
                row["requested_band"]
            ),
            source_band=np.asarray(
                row["source_band"]
            ),
            ra_deg=np.asarray(
                ra
            ),
            dec_deg=np.asarray(
                dec
            ),
        )


# ------------------------------------------------------------
# Build IRG / ZIR from the same per-band exposure products.
# ------------------------------------------------------------

def normalize_science_array(
    data,
    asinh_a=0.002,
    vmin=-0.03,
    vmax=3000,
):
    data = np.asarray(
        data,
        dtype=float,
    )

    norm = ImageNormalize(
        data,
        vmin=vmin,
        vmax=vmax,
        stretch=AsinhStretch(
            a=asinh_a
        ),
        clip=True,
    )

    return np.asarray(
        norm(data),
        dtype=float,
    )

def combine_rgb_arrays(
    red,
    green,
    blue,
):
    return np.dstack(
        [
            normalize_science_array(
                red
            ),
            normalize_science_array(
                green
            ),
            normalize_science_array(
                blue
            ),
        ]
    )


def rgb_array_to_png_bytes(
    rgb,
    saturation=RGB_SATURATION,
):
    arr = np.clip(
        rgb,
        0,
        1,
    )

    arr = np.uint8(
        np.round(
            255 * arr
        )
    )

    arr = np.flipud(
        arr
    )

    image = Image.fromarray(
        arr,
        mode="RGB",
    )

    if saturation != 1.0:
        enhancer = ImageEnhance.Color(
            image
        )
        image = enhancer.enhance(
            saturation
        )

    buffer = io.BytesIO()
    image.save(
        buffer,
        format="PNG",
    )
    return buffer.getvalue()

def source_image_for_requested_band(
    requested_band
):
    source_band = requested_to_source_band[
        requested_band
    ]

    return (
        band_products
        .get(
            source_band,
            {},
        )
        .get(
            "image"
        )
    )


irg_png = None
zir_png = None


irg_arrays = [
    source_image_for_requested_band(
        band
    )
    for band
    in [
        "i",
        "r",
        "g",
    ]
]

if all(
    x is not None
    for x
    in irg_arrays
):
    irg_png = rgb_array_to_png_bytes(
        combine_rgb_arrays(
            irg_arrays[0],
            irg_arrays[1],
            irg_arrays[2],
        )
    )


zir_arrays = [
    source_image_for_requested_band(
        band
    )
    for band
    in [
        "z",
        "i",
        "r",
    ]
]

if all(
    x is not None
    for x
    in zir_arrays
):
    zir_png = rgb_array_to_png_bytes(
        combine_rgb_arrays(
            zir_arrays[0],
            zir_arrays[1],
            zir_arrays[2],
        )
    )


print()
print(
    "IRG available:",
    irg_png is not None,
)

print(
    "ZIR available:",
    zir_png is not None,
)

if SAVE_DECONV_INPUT_NPZ:
    print(
        "Deconvolution input directory:",
        DECONV_INPUT_DIR,
    )


# 
# ## 6b. Placeholder multiband deconvolution with stacked `(6, 32, 32)` tensors
# 
# This section now prepares the multiband deconvolution inputs exactly as:
# 
# - `observed_object_images`: `(6, 32, 32)` in `ugrizy` order
# - `observed_psf_images`: `(6, 32, 32)` in `ugrizy` order
# 
# The placeholder deconvolution function returns:
# 
# - `predicted_object_images`: `(6, 32, 32)` in `ugrizy` order
# - `predicted_psf_images`: `(6, 32, 32)` in `ugrizy` order
# 
# You can insert your exact deconvolution model/function between the observed and
# predicted arrays.
# 
# The deconvolved IRG/ZIR color images are built from `predicted_object_images`.
# 
# 
# When predicted arrays are available, the notebook also saves one compressed file:
# 
# ```text
# target_<TARGET_ID>_deconvolution_arrays.npz
# ```
# 
# containing:
# 
# - `observed_object_images`
# - `observed_psf_images`
# - `predicted_object_images`
# - `predicted_psf_images`
# - `band_order`
# - `requested_to_source_band`
# - `ra_deg`
# - `dec_deg`
# 

# In[10]:


DECONV_STACK_BANDS = list(
    "ugrizy"
)

DECONV_CROP_SIZE = 32


def center_crop_or_pad(
    array,
    size=32,
    fill_value=0.0,
):
    """
    Return a centered square array of shape (size, size).

    If the input is larger, take the central crop.
    If the input is smaller, pad symmetrically.
    """
    array = np.asarray(
        array,
        dtype=float,
    )

    if array.ndim != 2:
        raise ValueError(
            "Input array must be 2D."
        )

    ny, nx = array.shape

    out = np.full(
        (
            size,
            size,
        ),
        fill_value,
        dtype=float,
    )

    # Input slice
    y0_in = max(
        0,
        (ny - size) // 2,
    )
    x0_in = max(
        0,
        (nx - size) // 2,
    )

    y1_in = min(
        ny,
        y0_in + size,
    )
    x1_in = min(
        nx,
        x0_in + size,
    )

    crop = array[
        y0_in:y1_in,
        x0_in:x1_in,
    ]

    cy, cx = crop.shape

    # Output slice
    y0_out = (size - cy) // 2
    x0_out = (size - cx) // 2

    out[
        y0_out:y0_out + cy,
        x0_out:x0_out + cx,
    ] = crop

    return out



def remove_single_batch_axis(
    array,
):
    """
    Remove the specific leading batch dimension if the shape is (1, 6, 32, 32).
    """
    array = np.asarray(
        array,
        dtype=float,
    )

    if array.shape == (
        1,
        6,
        DECONV_CROP_SIZE,
        DECONV_CROP_SIZE,
    ):
        array = array[0]

    return array


def deconvolve_multiband(
    observed_object_images,
    observed_psf_images,
    **kwargs,
):
    """
    Placeholder for multiband deconvolution.

    Parameters
    ----------
    observed_object_images : ndarray
        Shape (6, 32, 32), in ugrizy order.

    observed_psf_images : ndarray
        Shape (6, 32, 32), in ugrizy order.

    **kwargs
        Future deconvolution parameters.

    Returns
    -------
    predicted_object_images : ndarray
        Shape (6, 32, 32), in ugrizy order.

    predicted_psf_images : ndarray
        Shape (6, 32, 32), in ugrizy order.
    """
    observed_object_images = np.asarray(
        observed_object_images,
        dtype=float,
    )

    observed_psf_images = np.asarray(
        observed_psf_images,
        dtype=float,
    )

    expected_shape = (
        6,
        DECONV_CROP_SIZE,
        DECONV_CROP_SIZE,
    )

    if observed_object_images.shape != expected_shape:
        raise ValueError(
            "observed_object_images must have shape "
            f"{expected_shape}, got "
            f"{observed_object_images.shape}"
        )

    if observed_psf_images.shape != expected_shape:
        raise ValueError(
            "observed_psf_images must have shape "
            f"{expected_shape}, got "
            f"{observed_psf_images.shape}"
        )

    # --------------------------------------------------------
    # PLACEHOLDER
    #
    # Replace the next two lines with your exact deconvolution
    # function / model. Keep the outputs as:
    #
    #   predicted_object_images : (6,32,32)
    #   predicted_psf_images    : (6,32,32)
    #
    # in ugrizy order.
    # --------------------------------------------------------

    #predicted_object_images = observed_object_images.copy()
    #predicted_psf_images = observed_psf_images.copy()

    predicted_object_images, predicted_psf_images = my_enh.enhance(
        observed_object_images, 
        observed_psf_images,
        fix_psf=False,
    )

    return (
        predicted_object_images,
        predicted_psf_images,
    )



def normalize_for_rgb(
    data,
):
    data = np.asarray(
        data,
        dtype=float,
    )

    finite = data[
        np.isfinite(
            data
        )
    ]

    if len(finite) == 0:
        raise ValueError(
            "Input image has no finite pixels."
        )

    norm = ImageNormalize(
        data,
        stretch=LinearStretch(),
        clip=True,
    )

    return np.asarray(
        norm(data),
        dtype=float,
    )


def make_rgb_array(
    red,
    green,
    blue,
):
    red_n = normalize_for_rgb(
        red,
    )

    green_n = normalize_for_rgb(
        green,
    )

    blue_n = normalize_for_rgb(
        blue,
    )

    return np.clip(
        np.dstack(
            [
                red_n,
                green_n,
                blue_n,
            ]
        ),
        0,
        1,
    )


def make_deconvolved_color_images(
    predicted_object_images,
):
    """
    Build IRG and ZIR from predicted_object_images with shape (6,32,32)
    in ugrizy order.
    """
    predicted_object_images = np.asarray(
        predicted_object_images,
        dtype=float,
    )

    expected_shape = (
        6,
        DECONV_CROP_SIZE,
        DECONV_CROP_SIZE,
    )

    if predicted_object_images.shape != expected_shape:
        raise ValueError(
            "predicted_object_images must have shape "
            f"{expected_shape}, got "
            f"{predicted_object_images.shape}"
        )

    band_to_idx = {
        band: i
        for i, band in enumerate(DECONV_STACK_BANDS)
    }

    irg_rgb = make_rgb_array(
        red=predicted_object_images[
            band_to_idx["i"]
        ],
        green=predicted_object_images[
            band_to_idx["r"]
        ],
        blue=predicted_object_images[
            band_to_idx["g"]
        ],
    )

    zir_rgb = make_rgb_array(
        red=predicted_object_images[
            band_to_idx["z"]
        ],
        green=predicted_object_images[
            band_to_idx["i"]
        ],
        blue=predicted_object_images[
            band_to_idx["r"]
        ],
    )

    return {
        "irg_rgb":
            irg_rgb,
        "zir_rgb":
            zir_rgb,
        "irg_png":
            rgb_array_to_png_bytes(
                irg_rgb
            ),
        "zir_png":
            rgb_array_to_png_bytes(
                zir_rgb
            ),
    }


# ------------------------------------------------------------
# Build stacked observed arrays in ugrizy order.
# ------------------------------------------------------------

deconv_row_by_band = {
    row["requested_band"]:
        row
    for row
    in deconv_band_rows
}

observed_object_list = []
observed_psf_list = []

for band in DECONV_STACK_BANDS:

    if band not in deconv_row_by_band:
        raise KeyError(
            f"Missing requested band in deconv_band_rows: {band}"
        )

    row = deconv_row_by_band[
        band
    ]

    image = row.get(
        "image"
    )

    psf = row.get(
        "psf"
    )

    if image is None:
        raise ValueError(
            f"Missing science image for band {band}"
        )

    if psf is None:
        raise ValueError(
            f"Missing PSF image for band {band}"
        )

    image_32 = center_crop_or_pad(
        image,
        size=DECONV_CROP_SIZE,
        fill_value=0.0,
    )

    psf_32 = center_crop_or_pad(
        psf,
        size=DECONV_CROP_SIZE,
        fill_value=0.0,
    )

    observed_object_list.append(
        image_32
    )

    observed_psf_list.append(
        psf_32
    )


observed_object_images = np.stack(
    observed_object_list,
    axis=0,
)

observed_psf_images = np.stack(
    observed_psf_list,
    axis=0,
)

for i in range(
    observed_psf_images.shape[0]
):
    psf_sum = np.nansum(
        observed_psf_images[
            i
        ]
    )

    if (
        not np.isfinite(
            psf_sum
        )
        or psf_sum <= 0
    ):
        raise ValueError(
            "Cropped/padded PSF has a non-positive "
            f"sum for band {DECONV_STACK_BANDS[i]}."
        )

    observed_psf_images[
        i
    ] = (
        observed_psf_images[
            i
        ]
        / psf_sum
    )

for i, band in enumerate(
    DECONV_STACK_BANDS
):
    row = deconv_row_by_band[
        band
    ]

    row[
        "image_png_display"
    ] = grayscale_array_to_png_bytes(
        observed_object_images[
            i
        ]
    )

    row[
        "psf_png_display"
    ] = psf_array_to_png_bytes(
        observed_psf_images[
            i
        ]
    )

predicted_object_images = None
predicted_psf_images = None
deconv_irg_png = None
deconv_zir_png = None
predicted_object_pngs = {}

if RUN_DECONVOLUTION:

    (
        predicted_object_images,
        predicted_psf_images,
    ) = deconvolve_multiband(
        observed_object_images,
        observed_psf_images,
    )

    predicted_object_images = remove_single_batch_axis(
        predicted_object_images
    )

    predicted_psf_images = remove_single_batch_axis(
        predicted_psf_images
    )

    deconv_color_products = make_deconvolved_color_images(
        predicted_object_images
    )

    deconv_irg_png = deconv_color_products[
        "irg_png"
    ]

    deconv_zir_png = deconv_color_products[
        "zir_png"
    ]

    predicted_object_pngs = {
        band: grayscale_array_to_png_bytes(
            predicted_object_images[i]
        )
        for i, band in enumerate(
            DECONV_STACK_BANDS
        )
    }

    print(
        "Placeholder multiband deconvolution ran successfully."
    )

else:

    print(
        "RUN_DECONVOLUTION=False, so deconvolved IRG/ZIR were not produced."
    )


print(
    "observed_object_images shape:",
    observed_object_images.shape,
)

print(
    "observed_psf_images shape:",
    observed_psf_images.shape,
)

if predicted_object_images is not None:
    print(
        "predicted_object_images shape:",
        predicted_object_images.shape,
    )

if predicted_psf_images is not None:
    print(
        "predicted_psf_images shape:",
        predicted_psf_images.shape,
    )

print(
    "DECONV_STACK_BANDS =",
    DECONV_STACK_BANDS,
)


# ------------------------------------------------------------
# Save the four stacked deconvolution arrays in one NPZ file.
# ------------------------------------------------------------

DECONV_STACK_FILE = (
    DECONV_INPUT_DIR
    / (
        f"target_{target_id_safe_for_files}"
        "_deconvolution_arrays.npz"
    )
)

if (
    SAVE_DECONV_INPUT_NPZ
    and predicted_object_images is not None
    and predicted_psf_images is not None
):

    DECONV_INPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    np.savez_compressed(
        DECONV_STACK_FILE,

        observed_object_images=np.asarray(
            observed_object_images,
            dtype=float,
        ),

        observed_psf_images=np.asarray(
            observed_psf_images,
            dtype=float,
        ),

        predicted_object_images=np.asarray(
            predicted_object_images,
            dtype=float,
        ),

        predicted_psf_images=np.asarray(
            predicted_psf_images,
            dtype=float,
        ),

        band_order=np.asarray(
            DECONV_STACK_BANDS
        ),

        requested_to_source_band=np.asarray(
            [
                requested_to_source_band[
                    band
                ]
                for band
                in DECONV_STACK_BANDS
            ]
        ),

        ra_deg=np.asarray(
            ra
        ),

        dec_deg=np.asarray(
            dec
        ),
    )

    print(
        "Saved stacked deconvolution arrays:"
    )

    print(
        DECONV_STACK_FILE
    )

else:

    print(
        "Stacked deconvolution NPZ was not saved because "
        "predicted arrays are unavailable."
    )




# ## 7. Fetch Legacy Survey cutout

# In[11]:


legacy_url = (
    "https://www.legacysurvey.org/viewer/jpeg-cutout"
    f"?ra={ra:.8f}"
    f"&dec={dec:.8f}"
    f"&layer={LEGACY_LAYER}"
    f"&pixscale={LEGACY_PIXSCALE}"
    f"&size={LEGACY_SIZE_PIX}"
)

legacy_viewer_url = (
    "https://www.legacysurvey.org/viewer"
    f"?ra={ra:.8f}"
    f"&dec={dec:.8f}"
    f"&layer={LEGACY_LAYER}"
    "&zoom=16"
)

legacy_png_or_jpeg = None

try:
    response = requests.get(
        legacy_url,
        timeout=REQUEST_TIMEOUT,
    )

    response.raise_for_status()

    legacy_png_or_jpeg = response.content

    print(
        f"Legacy Survey cutout fetched: {len(response.content):,} bytes"
    )

except Exception as exc:
    print(
        "Legacy Survey cutout failed:",
        exc,
    )


# 
# ## 8. Prepare alerts for ANTARES thumbnail retrieval
# 
# The DP2 `Visit` table is not queried.
# 
# Alerts are sorted from newest to oldest. By default the notebook tries
# **every retained alert**. For a quick test, set:
# 
# ```python
# MAX_THUMBNAIL_ALERTS = 10
# ```
# 
# near the top of the notebook.
# 

# In[12]:


thumbnail_alerts = (
    target_alerts
    .drop_duplicates(
        subset=[
            "alert_id"
        ],
        keep="first",
    )
    .sort_values(
        "mjd",
        ascending=False,
    )
    .reset_index(
        drop=True
    )
)

if MAX_THUMBNAIL_ALERTS is not None:

    thumbnail_alerts = thumbnail_alerts.head(
        int(
            MAX_THUMBNAIL_ALERTS
        )
    ).copy()


print(
    "Retained target alerts:",
    len(
        target_alerts
    ),
)

print(
    "Unique alerts to try for thumbnails:",
    len(
        thumbnail_alerts
    ),
)

display(
    thumbnail_alerts[
        [
            "locus_id",
            "alert_id",
            "mjd",
            "band",
            "visit",
        ]
    ]
)


# 
# ## 9. ANTARES science / template / difference thumbnail retrieval
# 
# For each alert, retrieve the complete thumbnail packet with `get_thumbnails()`.
# 
# The helper accepts several plausible ANTARES response layouts, but the preferred
# layout is a mapping with the three keys:
# 
# ```text
# science
# template
# difference
# ```
# 
# Each alert is handled independently so one missing/bad thumbnail does not stop
# the rest of the target history.
# 

# In[13]:


THUMBNAIL_KINDS = [
    "science",
    "template",
    "difference",
]


def canonical_lsst_alert_id(
    alert_id,
):
    alert_id = str(
        alert_id
    ).strip()

    if alert_id.startswith(
        "lsst:"
    ):
        return alert_id

    return (
        "lsst:"
        + alert_id
    )


def normalize_thumbnail_label(
    value,
):
    if value is None:
        return ""

    return (
        str(
            value
        )
        .strip()
        .lower()
        .replace(
            "-",
            "_",
        )
        .replace(
            " ",
            "_",
        )
    )


THUMBNAIL_ALIASES = {
    "science": {
        "science",
        "visit",
        "science_image",
        "visit_image",
        "cutout_science",
        "cutoutscience",
    },
    "template": {
        "template",
        "template_image",
        "cutout_template",
        "cutouttemplate",
    },
    "difference": {
        "difference",
        "diff",
        "difference_image",
        "diff_image",
        "cutout_difference",
        "cutoutdifference",
    },
}


def item_label_values(
    item,
):
    if isinstance(
        item,
        Mapping,
    ):

        fields = [
            "type",
            "kind",
            "name",
            "thumbnail_type",
            "thumbnailType",
            "image_type",
            "imageType",
        ]

        return [
            normalize_thumbnail_label(
                item.get(
                    field
                )
            )
            for field
            in fields
        ]

    return []


def canonical_thumbnail_item(
    thumbnails,
    kind,
):
    aliases = THUMBNAIL_ALIASES[
        kind
    ]

    if isinstance(
        thumbnails,
        Mapping,
    ):

        # Most useful form:
        # {"science": ..., "template": ..., "difference": ...}
        for key, value in thumbnails.items():

            if (
                normalize_thumbnail_label(
                    key
                )
                in aliases
            ):
                return value

        # Common nested containers.
        for container_key in [
            "thumbnails",
            "data",
            "results",
            "cutouts",
        ]:

            if container_key in thumbnails:

                try:
                    return canonical_thumbnail_item(
                        thumbnails[
                            container_key
                        ],
                        kind=kind,
                    )

                except Exception:
                    pass

        # Single labeled item.
        labels = item_label_values(
            thumbnails
        )

        if any(
            label in aliases
            for label
            in labels
        ):
            return thumbnails


    if isinstance(
        thumbnails,
        (
            list,
            tuple,
        ),
    ):

        for item in thumbnails:

            if not isinstance(
                item,
                Mapping,
            ):
                continue

            labels = item_label_values(
                item
            )

            if any(
                label in aliases
                for label
                in labels
            ):
                return item

        for item in thumbnails:

            if not isinstance(
                item,
                Mapping,
            ):
                continue

            for key, value in item.items():

                if (
                    normalize_thumbnail_label(
                        key
                    )
                    in aliases
                ):
                    return value


    raise KeyError(
        f"Could not identify ANTARES {kind!r} thumbnail."
    )


def thumbnail_item_to_png_bytes(
    item,
):
    """
    Convert one ANTARES thumbnail item to ordinary PNG bytes.

    Supports:
      - mapping with blob/file_bytes/bytes/data
      - direct bytes-like object
      - base64/data-URI string
      - URL string
      - PIL Image
      - numpy image array
    """

    if item is None:
        return None


    if isinstance(
        item,
        Image.Image,
    ):

        image = item.convert(
            "RGB"
        )

        buffer = io.BytesIO()

        image.save(
            buffer,
            format="PNG",
        )

        return buffer.getvalue()


    if isinstance(
        item,
        np.ndarray,
    ):

        data = np.asarray(
            item
        )

        if data.dtype != np.uint8:

            finite = data[
                np.isfinite(
                    data
                )
            ]

            if len(
                finite
            ) == 0:
                return None

            lo, hi = np.nanpercentile(
                finite,
                [
                    1,
                    99,
                ],
            )

            if hi <= lo:
                return None

            data = np.clip(
                (
                    data
                    - lo
                )
                / (
                    hi
                    - lo
                ),
                0,
                1,
            )

            data = np.uint8(
                np.round(
                    255
                    * np.nan_to_num(
                        data,
                        nan=0.0,
                    )
                )
            )

        image = Image.fromarray(
            data
        ).convert(
            "RGB"
        )

        buffer = io.BytesIO()

        image.save(
            buffer,
            format="PNG",
        )

        return buffer.getvalue()


    value = item

    if isinstance(
        item,
        Mapping,
    ):

        value = None

        for key in [
            "blob",
            "file_bytes",
            "bytes",
            "data",
            "image",
        ]:

            if key in item:

                value = item[
                    key
                ]

                break

        if value is None:

            raise KeyError(
                "Thumbnail item has no recognized image payload field."
            )


    if isinstance(
        value,
        memoryview,
    ):
        value = value.tobytes()

    elif isinstance(
        value,
        bytearray,
    ):
        value = bytes(
            value
        )


    if isinstance(
        value,
        str,
    ):

        text = value.strip()

        if (
            text.startswith(
                "http://"
            )
            or text.startswith(
                "https://"
            )
        ):

            response = requests.get(
                text,
                timeout=REQUEST_TIMEOUT,
            )

            response.raise_for_status()

            value = response.content

        else:

            if text.startswith(
                "data:"
            ):
                text = text.split(
                    ",",
                    1,
                )[1]

            value = base64.b64decode(
                text
            )


    if not isinstance(
        value,
        bytes,
    ):

        raise TypeError(
            "Thumbnail payload is not bytes-like: "
            f"{type(value)}"
        )


    image = Image.open(
        io.BytesIO(
            value
        )
    )

    image.load()

    image = image.convert(
        "RGB"
    )

    buffer = io.BytesIO()

    image.save(
        buffer,
        format="PNG",
    )

    return buffer.getvalue()


def fetch_antares_thumbnail_triplet(
    alert_id,
):
    """
    Retrieve science, template, and difference thumbnails for one alert.

    Returns
    -------
    products : dict
        PNG bytes or None for each of:
        science, template, difference

    errors : dict
        Per-kind decode/lookup failures.
    """

    antares_alert_id = canonical_lsst_alert_id(
        alert_id
    )

    thumbnails = get_thumbnails(
        antares_alert_id
    )

    products = {}
    errors = {}


    for kind in THUMBNAIL_KINDS:

        try:

            item = canonical_thumbnail_item(
                thumbnails,
                kind=kind,
            )

            products[
                kind
            ] = thumbnail_item_to_png_bytes(
                item
            )

        except Exception as exc:

            products[
                kind
            ] = None

            errors[
                kind
            ] = str(
                exc
            )


    return (
        products,
        errors,
        thumbnails,
    )


# 
# ## 10. Fetch as many ANTARES alert thumbnail triplets as possible
# 
# The notebook attempts every selected alert, newest first.
# 
# For each alert it keeps:
# 
# - science / visit thumbnail
# - template thumbnail
# - difference thumbnail
# - alert ID
# - MJD
# - band
# - visit number
# - locus ID
# 
# No thumbnail is rotated.
# 

# In[14]:


thumbnail_frames = []
thumbnail_failures = []


n_to_try = len(
    thumbnail_alerts
)


for i, row in thumbnail_alerts.iterrows():

    locus_id = str(
        row[
            "locus_id"
        ]
    )

    alert_id = str(
        row[
            "alert_id"
        ]
    )

    try:

        products, per_kind_errors, raw_thumbnails = fetch_antares_thumbnail_triplet(
            alert_id
        )

        n_available = sum(
            products[
                kind
            ]
            is not None
            for kind
            in THUMBNAIL_KINDS
        )

        if THUMBNAIL_VERBOSE:

            print(
                f"[{i + 1}/{n_to_try}] {alert_id}: "
                f"science={products['science'] is not None}, "
                f"template={products['template'] is not None}, "
                f"difference={products['difference'] is not None}"
            )


        if n_available == 0:

            thumbnail_failures.append(
                {
                    "alert_id":
                        alert_id,
                    "error":
                        (
                            "No science/template/difference thumbnails "
                            "could be decoded."
                        ),
                    "per_kind_errors":
                        per_kind_errors,
                }
            )

            continue


        thumbnail_frames.append(
            {
                "locus_id":
                    locus_id,

                "alert_id":
                    alert_id,

                "mjd":
                    float(
                        row[
                            "mjd"
                        ]
                    ),

                "band":
                    (
                        None
                        if pd.isna(
                            row[
                                "band"
                            ]
                        )
                        else str(
                            row[
                                "band"
                            ]
                        )
                    ),

                "visit":
                    row[
                        "visit"
                    ],

                "science_png":
                    products[
                        "science"
                    ],

                "template_png":
                    products[
                        "template"
                    ],

                "difference_png":
                    products[
                        "difference"
                    ],

                "per_kind_errors":
                    per_kind_errors,
            }
        )


    except Exception as exc:

        if THUMBNAIL_VERBOSE:

            print(
                f"[{i + 1}/{n_to_try}] {alert_id}: FAILED: {exc}"
            )

        thumbnail_failures.append(
            {
                "alert_id":
                    alert_id,
                "error":
                    str(
                        exc
                    ),
                "per_kind_errors":
                    {},
            }
        )


n_complete_triplets = sum(
    (
        frame[
            "science_png"
        ]
        is not None

        and frame[
            "template_png"
        ]
        is not None

        and frame[
            "difference_png"
        ]
        is not None
    )
    for frame
    in thumbnail_frames
)


print()
print(
    "Thumbnail retrieval summary:"
)

print(
    f"  tried alerts:              {n_to_try}"
)

print(
    f"  alerts with >=1 thumbnail: {len(thumbnail_frames)}"
)

print(
    f"  complete triplets:         {n_complete_triplets}"
)

print(
    f"  failed/no-image alerts:    {len(thumbnail_failures)}"
)


# ## 11. HTML helpers

# In[15]:


def fmt(
    value,
    digits=3,
):
    try:
        if pd.isna(
            value
        ):
            return "N/A"
    except Exception:
        pass

    try:
        return (
            f"{float(value):.{digits}f}"
        )
    except Exception:
        return html.escape(
            str(
                value
            )
        )


def bool_text(
    value
):
    try:
        if pd.isna(
            value
        ):
            return "N/A"
    except Exception:
        pass

    return (
        "Yes"
        if bool(
            value
        )
        else "No"
    )


def bytes_to_data_uri(
    image_bytes,
    mime="image/png",
):
    if image_bytes is None:
        return None

    encoded = base64.b64encode(
        image_bytes
    ).decode(
        "ascii"
    )

    return (
        f"data:{mime};base64,{encoded}"
    )


def placeholder_card(
    title,
    message="Not available",
):
    return f"""
    <div class="image-card">
      <div class="image-title">{html.escape(title)}</div>
      <div class="placeholder">{html.escape(message)}</div>
    </div>
    """


def image_card(
    title,
    image_uri,
    caption="",
    link=None,
    pixelated=False,
):
    if image_uri is None:
        return placeholder_card(
            title
        )

    image_class = (
        "pixelated"
        if pixelated
        else ""
    )

    image_html = (
        f'<img class="{image_class}" '
        f'src="{image_uri}" '
        f'alt="{html.escape(title)}">'
    )

    if link:
        image_html = (
            f'<a href="{html.escape(link)}" target="_blank">'
            f'{image_html}'
            f'</a>'
        )

    caption_html = (
        f'<div class="caption">{caption}</div>'
        if caption
        else ""
    )

    return f"""
    <div class="image-card">
      <div class="image-title">{html.escape(title)}</div>
      {image_html}
      {caption_html}
    </div>
    """


def alert_triplet_html(
    frame,
):
    metadata = (
        f"MJD {frame['mjd']:.5f}"
        f" | band {html.escape(str(frame['band']))}"
        f" | visit {html.escape(str(frame['visit']))}"
        f" | alert {html.escape(frame['alert_id'])}"
        f" | locus {html.escape(frame['locus_id'])}"
    )

    science_uri = bytes_to_data_uri(
        frame[
            "science_png"
        ]
    )

    template_uri = bytes_to_data_uri(
        frame[
            "template_png"
        ]
    )

    difference_uri = bytes_to_data_uri(
        frame[
            "difference_png"
        ]
    )

    cards = "\n".join(
        [
            image_card(
                "Science / Visit",
                science_uri,
                pixelated=True,
            ),
            image_card(
                "Template",
                template_uri,
                pixelated=True,
            ),
            image_card(
                "Difference",
                difference_uri,
                pixelated=True,
            ),
        ]
    )

    return f"""
    <div class="alert-block">
      <div class="alert-metadata">{metadata}</div>
      <div class="image-grid alert-triplet">
        {cards}
      </div>
    </div>
    """

def deconv_row_html(
    requested_band,
    source_band,
    image_uri,
    psf_uri,
    deconvolved_uri,
    used_fallback=None,
):
    if used_fallback is None:
        used_fallback = (
            requested_band
            != source_band
        )

    fallback_text = (
        " [fallback]"
        if used_fallback
        else ""
    )

    if used_fallback:
        band_label = (
            f"{requested_band} → {source_band} [fallback]"
        )
    else:
        band_label = (
            f"{requested_band} band"
        )

    science_card = image_card(
        f"{band_label}: science cutout",
        image_uri,
        pixelated=True,
    )

    psf_card = image_card(
        f"{band_label}: PSF kernel",
        psf_uri,
        pixelated=True,
    )

    deconv_card = image_card(
        f"{band_label}: deconvolved",
        deconvolved_uri,
        pixelated=True,
    )

    return f"""
    <div class="deconv-band-title">{html.escape(band_label)}</div>
    <div class="image-grid deconv-row">
      {science_card}
      {psf_card}
      {deconv_card}
    </div>
    """


# ## 12. Build the self-contained HTML page

# In[16]:


# ------------------------------------------------------------
# Flags and scalar metadata
# ------------------------------------------------------------

gaia_qso = target.get(
    "gaia_in_qso_candidates_any",
    np.nan,
)

milliquas = target.get(
    "milliquas_match",
    np.nan,
)

gaia_bp_rp = target.get(
    "gaia_bp_rp_max",
    np.nan,
)

gaia_pm_snr = target.get(
    "gaia_pm_sig_max",
    np.nan,
)

gaia_parallax_snr = target.get(
    "gaia_parallax_sig_max",
    np.nan,
)

second_pass_score = target.get(
    "second_pass_score",
    np.nan,
)

first_cut_score = target.get(
    "first_cut_max_score",
    np.nan,
)


if known_lens_match is None:
    known_lens_text = "N/A"

elif known_lens_match:

    known_lens_text = "Yes"

    if known_lens_name:

        known_lens_text += (
            f" ({html.escape(known_lens_name)})"
        )

else:
    known_lens_text = "No"


# ------------------------------------------------------------
# Magnitude HTML
# ------------------------------------------------------------

mag_cells = []

for band in "ugrizy":

    mag_cells.append(
        f"""
        <div class="metric">
          <div class="metric-label">{band}</div>
          <div class="metric-value">{fmt(latest_alert_mags[band], 2)}</div>
        </div>
        """
    )

magnitudes_html = "\n".join(
    mag_cells
)


# ------------------------------------------------------------
# Band counts
# ------------------------------------------------------------

band_count_parts = []

for band in "ugrizy":

    n = int(
        band_counts.get(
            band,
            0,
        )
    )

    if n > 0:

        band_count_parts.append(
            f"{band}: {n}"
        )


band_count_text = (
    ", ".join(
        band_count_parts
    )
    if band_count_parts
    else "None"
)


# ------------------------------------------------------------
# Static image cards
# ------------------------------------------------------------

irg_uri = bytes_to_data_uri(
    irg_png
)

zir_uri = bytes_to_data_uri(
    zir_png
)

legacy_uri = bytes_to_data_uri(
    legacy_png_or_jpeg,
    mime="image/jpeg",
)


static_cards = "\n".join(
    [
        image_card(
            "IRG",
            irg_uri,
            pixelated=True,
        ),

        image_card(
            "ZIR",
            zir_uri,
            pixelated=True,
        ),

        image_card(
            "Legacy Survey",
            legacy_uri,
            caption=(
                '<a href="'
                + html.escape(
                    legacy_viewer_url
                )
                + '" target="_blank">'
                + "Open Legacy Survey viewer"
                + "</a>"
            ),
            link=legacy_viewer_url,
            pixelated=True,
        ),
    ]
)


# ------------------------------------------------------------
# Deconvolution-preparation rows
# ------------------------------------------------------------

deconv_rows_html = []

for row in deconv_band_rows:

    requested_band = row[
        "requested_band"
    ]

    image_uri = bytes_to_data_uri(
        row.get(
            "image_png_display",
            row[
                "image_png"
            ],
        )
    )

    psf_uri = bytes_to_data_uri(
        row.get(
            "psf_png_display",
            row[
                "psf_png"
            ],
        )
    )

    deconvolved_uri = bytes_to_data_uri(
        predicted_object_pngs.get(
            requested_band
        )
    )

    deconv_rows_html.append(
        deconv_row_html(
            requested_band=requested_band,
            source_band=row[
                "source_band"
            ],
            image_uri=image_uri,
            psf_uri=psf_uri,
            deconvolved_uri=deconvolved_uri,
            used_fallback=(
                requested_band
                != row[
                    "source_band"
                ]
            ),
        )
    )


deconv_rows_html = "\n".join(
    deconv_rows_html
)


band_mapping_text = ", ".join(
    (
        f"{row['requested_band']}→{row['source_band']}"
        if row[
            "requested_band"
        ]
        != row[
            "source_band"
        ]
        else row[
            "requested_band"
        ]
    )
    for row
    in deconv_band_rows
)


fallback_rows = [
    row
    for row
    in deconv_band_rows
    if row[
        "requested_band"
    ] != row[
        "source_band"
    ]
]

if fallback_rows:

    fallback_summary_text = ", ".join(
        f"{row['requested_band']}→{row['source_band']}"
        for row
        in fallback_rows
    )

    fallback_html = (
        "<div class=\"metric\">"
        "<div class=\"metric-label\">Fallback replaced bands</div>"
        f"<div class=\"metric-value\">{html.escape(fallback_summary_text)}</div>"
        "</div>"
    )

    fallback_note_html = (
        "<div class=\"subtle\" style=\"margin-bottom:10px;\">"
        "Some requested bands were unavailable and were replaced by the nearest usable band. "
        "Rows marked with <code>[fallback]</code> indicate those replacements."
        "</div>"
    )

else:

    fallback_summary_text = "None"

    fallback_html = (
        "<div class=\"metric\">"
        "<div class=\"metric-label\">Fallback replaced bands</div>"
        "<div class=\"metric-value\">None</div>"
        "</div>"
    )

    fallback_note_html = (
        "<div class=\"subtle\" style=\"margin-bottom:10px;\">"
        "No requested bands needed fallback replacement."
        "</div>"
    )


# ------------------------------------------------------------
# Deconvolved color composites
# ------------------------------------------------------------

deconv_irg_uri = bytes_to_data_uri(
    deconv_irg_png
)

deconv_zir_uri = bytes_to_data_uri(
    deconv_zir_png
)


if (
    RUN_DECONVOLUTION
    and deconv_irg_uri is not None
    and deconv_zir_uri is not None
):

    deconv_color_cards = "\n".join(
        [
            image_card(
                "Deconvolved IRG",
                deconv_irg_uri,
                pixelated=True,
            ),

            image_card(
                "Deconvolved ZIR",
                deconv_zir_uri,
                pixelated=True,
            ),
        ]
    )

else:

    deconv_color_cards = "\n".join(
        [
            placeholder_card(
                "Deconvolved IRG",
                "Not produced yet.",
            ),
            placeholder_card(
                "Deconvolved ZIR",
                "Not produced yet.",
            ),
        ]
    )


# ------------------------------------------------------------
# All successful/partial thumbnail triplets
# ------------------------------------------------------------

alert_triplets_html = "\n".join(
    alert_triplet_html(
        frame
    )
    for frame
    in thumbnail_frames
)


# ------------------------------------------------------------
# JavaScript animation payload.
# Keep science/template/difference synchronized by alert.
# ------------------------------------------------------------

animation_frames = []

for frame in thumbnail_frames:

    animation_frames.append(
        {
            "science":
                bytes_to_data_uri(
                    frame[
                        "science_png"
                    ]
                ),

            "template":
                bytes_to_data_uri(
                    frame[
                        "template_png"
                    ]
                ),

            "difference":
                bytes_to_data_uri(
                    frame[
                        "difference_png"
                    ]
                ),

            "caption":
                (
                    f"MJD {frame['mjd']:.5f}"
                    f" | band {frame['band']}"
                    f" | visit {frame['visit']}"
                    f" | alert {frame['alert_id']}"
                ),
        }
    )


animation_frames_json = json.dumps(
    animation_frames
)


# ------------------------------------------------------------
# Write HTML
# ------------------------------------------------------------

target_id_safe = re.sub(
    r"[^A-Za-z0-9_.-]+",
    "_",
    str(
        TARGET_ID
    ),
)

html_path = (
    OUTPUT_DIR
    / f"{target_id_safe}.html"
)

target_summary_line = ",".join(
    [
        str(TARGET_ID),
        f"{ra:.8f}",
        f"{dec:.8f}",
    ]
    + [
        fmt(
            latest_alert_mags[band],
            4,
        )
        for band in "ugrizy"
    ]
)

page = f"""<!doctype html>
<html lang="en">

<head>

<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">

<title>{html.escape(str(TARGET_ID))}</title>

<style>

  :root {{
    color-scheme: light dark;
    --bg: #111318;
    --panel: #1b1f27;
    --panel2: #242a34;
    --text: #f2f4f7;
    --muted: #aeb7c4;
    --border: #394250;
    --accent: #7cb8ff;
  }}

  body {{
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    background: var(--bg);
    color: var(--text);
  }}

  .page {{
    max-width: 1500px;
    margin: 0 auto;
    padding: 24px;
  }}

  h1 {{
    margin: 0 0 8px 0;
    font-size: 28px;
  }}

  h2 {{
    margin-top: 32px;
    font-size: 20px;
  }}

  .subtle {{
    color: var(--muted);
  }}

  .summary-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
    gap: 10px;
    margin-top: 16px;
  }}

  .metric {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 10px 12px;
  }}

  .metric-label {{
    color: var(--muted);
    font-size: 12px;
  }}

  .metric-value {{
    margin-top: 3px;
    font-size: 16px;
    word-break: break-word;
  }}

  .mag-grid {{
    display: grid;
    grid-template-columns: repeat(6, minmax(80px, 1fr));
    gap: 8px;
  }}

  .image-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
    gap: 14px;
  }}

  .image-grid.two {{
    grid-template-columns: repeat(2, minmax(280px, 1fr));
  }}

  .deconv-row,
  .alert-triplet,
  .animation-triplet {{
    grid-template-columns: repeat(3, minmax(240px, 1fr));
  }}

  .image-card {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 10px;
    text-align: center;
  }}

  .image-title {{
    font-weight: 650;
    margin-bottom: 8px;
  }}

  .image-card img,
  .animation-card img {{
    display: block;
    width: 100%;
    max-height: 520px;
    object-fit: contain;
    background: black;
    border-radius: 6px;
  }}

  img.pixelated,
  .alert-triplet img,
  .animation-triplet img {{
    image-rendering: pixelated;
    image-rendering: crisp-edges;
  }}

  .caption {{
    margin-top: 8px;
    color: var(--muted);
    font-size: 12px;
    line-height: 1.45;
  }}

  .placeholder {{
    min-height: 240px;
    display: flex;
    align-items: center;
    justify-content: center;
    border: 1px dashed var(--border);
    border-radius: 6px;
    color: var(--muted);
    padding: 12px;
  }}

  .deconv-band-title {{
    margin-top: 18px;
    margin-bottom: 8px;
    font-weight: 700;
    font-size: 16px;
  }}

  .alert-block {{
    margin-top: 18px;
    padding-top: 14px;
    border-top: 1px solid var(--border);
  }}

  .alert-metadata {{
    margin-bottom: 8px;
    font-weight: 650;
    font-size: 14px;
  }}

  .animation {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 14px;
  }}

  .animation-card {{
    text-align: center;
  }}

  .animation-title {{
    font-weight: 650;
    margin-bottom: 8px;
  }}

  .animation-missing {{
    min-height: 260px;
    display: none;
    align-items: center;
    justify-content: center;
    border: 1px dashed var(--border);
    border-radius: 6px;
    color: var(--muted);
  }}

  .controls {{
    margin-top: 12px;
    display: flex;
    gap: 8px;
    justify-content: center;
    flex-wrap: wrap;
  }}

  button {{
    padding: 7px 14px;
    border: 1px solid var(--border);
    border-radius: 7px;
    background: var(--panel2);
    color: var(--text);
    cursor: pointer;
  }}

  a {{
    color: var(--accent);
  }}

  @media (max-width: 800px) {{

    .deconv-row,
    .alert-triplet,
    .animation-triplet {{
      grid-template-columns: 1fr;
    }}

  }}

</style>

</head>


<body>

<div class="page">


<h1>{html.escape(str(TARGET_ID))}</h1>

<div class="subtle" style="margin-bottom:12px; font-family:monospace;">
  {html.escape(target_summary_line)}
</div>


<h2>Static imaging</h2>

<div class="image-grid">
{static_cards}
</div>


<div class="summary-grid">

  <div class="metric">
    <div class="metric-label">RA</div>
    <div class="metric-value">{ra:.8f}°</div>
  </div>

  <div class="metric">
    <div class="metric-label">Dec</div>
    <div class="metric-value">{dec:.8f}°</div>
  </div>

  <div class="metric">
    <div class="metric-label">Galactic latitude b</div>
    <div class="metric-value">{galactic_latitude_deg:.4f}°</div>
  </div>

  <div class="metric">
    <div class="metric-label">Latest-alert locus ID</div>
    <div class="metric-value">{html.escape(latest_locus_id)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">Second-Pass score</div>
    <div class="metric-value">{fmt(second_pass_score, 4)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">First-Cut score</div>
    <div class="metric-value">{fmt(first_cut_score, 4)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">Gaia PM S/N (max)</div>
    <div class="metric-value">{fmt(gaia_pm_snr, 2)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">Gaia parallax S/N (max)</div>
    <div class="metric-value">{fmt(gaia_parallax_snr, 2)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">Gaia BP-RP (max)</div>
    <div class="metric-value">{fmt(gaia_bp_rp, 2)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">GaiaQSO</div>
    <div class="metric-value">{bool_text(gaia_qso)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">Milliquas</div>
    <div class="metric-value">{bool_text(milliquas)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">KnownLens</div>
    <div class="metric-value">{known_lens_text}</div>
  </div>

  <div class="metric">
    <div class="metric-label">DiaSource extendedness median</div>
    <div class="metric-value">{fmt(extendedness_median, 3)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">N alert after 2026-05-27</div>
    <div class="metric-value">{len(target_alerts)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">Bands / counts</div>
    <div class="metric-value">{html.escape(band_count_text)}</div>
  </div>

  <div class="metric">
    <div class="metric-label">Thumbnail alerts retrieved</div>
    <div class="metric-value">{len(thumbnail_frames)} / {len(thumbnail_alerts)}</div>
  </div>

  {fallback_html}

</div>


<h2>Most-recent-alert diaObject magnitudes</h2>

<div class="subtle" style="margin-bottom:10px;">
  AB magnitudes from lsst_diaObject_&#123;band&#125;_scienceFluxMean
  in the most recent retained alert
  (MJD {float(latest_alert['mjd']):.5f},
  {html.escape(str(latest_alert['band']))} band).
</div>

<div class="mag-grid">
{magnitudes_html}
</div>



<h2>Deconvolved color composites</h2>

<div class="image-grid two">
{deconv_color_cards}
</div>


<h2>Deconvolution preparation by band</h2>

<div class="subtle" style="margin-bottom:10px;">
  Band mapping (if unavailable): {html.escape(band_mapping_text)}.
  The science cutout and PSF kernel come from the DP2 deep coadd.
  Science cutout, PSF kernel, and deconvolved images are displayed
  using a linear stretch.
</div>

{deconv_rows_html}


<h2>Animated ANTARES alert thumbnails</h2>

<div class="animation">

  <div class="image-grid animation-triplet">

    <div class="animation-card">
      <div class="animation-title">Science / Visit</div>
      <img id="animScience" class="pixelated" alt="science thumbnail">
      <div id="missingScience" class="animation-missing">Not available</div>
    </div>

    <div class="animation-card">
      <div class="animation-title">Template</div>
      <img id="animTemplate" class="pixelated" alt="template thumbnail">
      <div id="missingTemplate" class="animation-missing">Not available</div>
    </div>

    <div class="animation-card">
      <div class="animation-title">Difference</div>
      <img id="animDifference" class="pixelated" alt="difference thumbnail">
      <div id="missingDifference" class="animation-missing">Not available</div>
    </div>

  </div>

  <div id="animCaption" class="caption"></div>

  <div class="controls">
    <button onclick="previousFrame()">Previous</button>
    <button id="playButton" onclick="togglePlay()">Pause</button>
    <button onclick="nextFrame()">Next</button>
  </div>

</div>


<h2>All retrieved ANTARES alert thumbnails</h2>

<div class="subtle" style="margin-bottom:10px;">
  Alerts are ordered newest to oldest.
</div>

{alert_triplets_html}


</div>


<script>

const frames = {animation_frames_json};
const intervalMs = {int(ANIMATION_INTERVAL_MS)};

let frameIndex = 0;
let playing = true;
let timer = null;


function setImageOrMissing(imageId, missingId, src) {{

  const image =
    document.getElementById(imageId);

  const missing =
    document.getElementById(missingId);

  if (src) {{

    image.src = src;
    image.style.display = "block";
    missing.style.display = "none";

  }}

  else {{

    image.removeAttribute("src");
    image.style.display = "none";
    missing.style.display = "flex";

  }}

}}


function showFrame(index) {{

  if (frames.length === 0) {{

    document.getElementById("animCaption").textContent =
      "No ANTARES alert thumbnails were retrieved.";

    setImageOrMissing(
      "animScience",
      "missingScience",
      null
    );

    setImageOrMissing(
      "animTemplate",
      "missingTemplate",
      null
    );

    setImageOrMissing(
      "animDifference",
      "missingDifference",
      null
    );

    return;

  }}


  frameIndex =
    (index + frames.length)
    % frames.length;


  const frame =
    frames[
      frameIndex
    ];


  setImageOrMissing(
    "animScience",
    "missingScience",
    frame.science
  );

  setImageOrMissing(
    "animTemplate",
    "missingTemplate",
    frame.template
  );

  setImageOrMissing(
    "animDifference",
    "missingDifference",
    frame.difference
  );


  document.getElementById("animCaption").textContent =
    `${{frameIndex + 1}} / ${{frames.length}} · ${{frame.caption}}`;

}}


function nextFrame() {{

  showFrame(
    frameIndex + 1
  );

}}


function previousFrame() {{

  showFrame(
    frameIndex - 1
  );

}}


function stopTimer() {{

  if (timer !== null) {{

    clearInterval(
      timer
    );

    timer = null;

  }}

}}


function startTimer() {{

  stopTimer();

  if (
    playing
    && frames.length > 1
  ) {{

    timer = setInterval(
      nextFrame,
      intervalMs
    );

  }}

}}


function togglePlay() {{

  playing =
    !playing;


  const button =
    document.getElementById(
      "playButton"
    );


  if (playing) {{

    button.textContent =
      "Pause";

    startTimer();

  }}

  else {{

    button.textContent =
      "Play";

    stopTimer();

  }}

}}


showFrame(
  0
);

startTimer();

</script>


</body>

</html>
"""


html_path.write_text(
    page,
    encoding="utf-8",
)


print(
    "HTML written:"
)

print(
    html_path
)

display(
    FileLink(
        str(
            html_path
        )
    )
)


# 
# ## v21 test checklist
# 
# 1. **Top settings**
#    - run the first settings cell and enter `TARGET_ID` when prompted;
#    - leave `MAX_THUMBNAIL_ALERTS = None` to try all alerts;
#    - for a quick test, use an integer such as `10`;
# 
# 2. **Target metadata**
#    - RA/Dec and Galactic latitude;
#    - latest-alert locus ID;
#    - First-Cut / Second-Pass scores;
#    - Gaia PM/parallax significance;
#    - GaiaQSO, Milliquas, KnownLens.
# 
# 3. **DiaSource extendedness**
#    - confirm the median is computed from retained `lsst_diaSource_extendedness` values;
#    - the value is computed only from
#      `lsst_diaSource_extendedness` in retained alerts.
# 
# 4. **Static imaging**
#    - IRG / ZIR should show native pixels rather than browser-smoothed pixels;
#    - adjust `RGB_SATURATION` near the top of the notebook to change fake-color saturation;
#    - Legacy Survey is displayed with nearest-neighbor / pixelated rendering.
# 
# 5. **Deconvolution tensors**
#    - confirm
#      `observed_object_images.shape == (6, 32, 32)`;
#    - confirm
#      `observed_psf_images.shape == (6, 32, 32)`;
#    - confirm, if `RUN_DECONVOLUTION=True`,
#      `predicted_object_images.shape == (6, 32, 32)` and
#      `predicted_psf_images.shape == (6, 32, 32)`;
#    - confirm the band order is `ugrizy`;
#    - insert the real deconvolution function inside
#      `deconvolve_multiband(...)`;
#    - confirm `target_<TARGET_ID>_deconvolution_arrays.npz` is saved
#      and contains all four `(6, 32, 32)` arrays.
# 
# 6. **Deconvolution preparation / missing bands**
#    - confirm `usable_bands` contains only bands whose actual cutout and PSF succeeded;
#    - confirm the final fallback mapping is based on `usable_bands`, not only the SIA listing;
#    - for example, if `u` fails but `g` works, confirm `u -> g [fallback]`;
#    - confirm the HTML summary shows the replaced bands;
#    - confirm the per-band rows show `[fallback]` where replacement happened;
#    - inspect science cutouts and PSF kernels;
#    - confirm saved `.npz` inputs;
#    - confirm the per-band `deconvolved` panels are filled from `predicted_object_images`.
# 
# 7. **ANTARES alert thumbnails**
#    - no `Visit` table is queried;
#    - no rotation is applied;
#    - science/visit, template, and difference are shown together;
#    - the notebook tries all retained alerts unless capped;
#    - retrieval failures should not stop later alerts.
# 
# 8. **Animation**
#    - science/template/difference advance synchronously;
#    - play / pause / previous / next work;
#    - alert metadata changes with each frame.
# 
