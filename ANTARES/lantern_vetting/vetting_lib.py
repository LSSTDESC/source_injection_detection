"""
Lens-vetting helper library for ANTARES loci.

Cutout fetchers, catalog queries, scale bar, ANTARES search,
Gaia stellar filter, and CSV classification I/O.
"""

import csv
import os
from datetime import datetime
from io import BytesIO
from itertools import islice

import numpy as np
import requests
from PIL import Image, ImageDraw
from antares_client.search import search, get_by_id, get_thumbnails
from astroquery.mast import Observations
from astroquery.esa.euclid import Euclid
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u


# --- Cutout fetchers -------------------------------------------------------

def get_decals_jpg(ra, dec, size=100, layer="ls-dr10-grz", pixscale=0.262):
    """Fetch DECaLS DR10 colour JPG cutout."""
    url = (f"https://www.legacysurvey.org/viewer/cutout.jpg"
           f"?ra={ra}&dec={dec}&layer={layer}&pixscale={pixscale}&size={size}")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.content


def get_ps1_color_cutout(ra, dec, size=120, output_size=256):
    """Fetch PS1 color cutout. Prefer grz, fall back to gri.
    size = cutout in PS1 pixels (0.25"/pix), so 120 = 30 arcsec on sky.
    """
    url = (f"https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
           f"?ra={ra}&dec={dec}&filters=griz&type=stack")
    lines = requests.get(url, timeout=30).text.strip().split('\n')
    files = {}
    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 8:
            files[parts[4]] = parts[7]
    if all(f in files for f in 'grz'):
        red, green, blue, label = files['z'], files['r'], files['g'], 'grz'
    elif all(f in files for f in 'gri'):
        red, green, blue, label = files['i'], files['r'], files['g'], 'gri'
    else:
        return None, None
    url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
           f"?ra={ra}&dec={dec}&size={size}&format=jpg&output_size={output_size}"
           f"&red={red}&green={green}&blue={blue}")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.content, label


def get_hst_cutout(ra, dec, fov_deg=0.003, size=256):
    """Fetch HST cutout via CDS HiPS2FITS, trying multiple band maps."""
    for hips in ['CDS/P/HST/wideV', 'CDS/P/HST/color', 'CDS/P/HST/I',
                 'CDS/P/HST/R', 'CDS/P/HST/V', 'CDS/P/HST/B',
                 'CDS/P/HST/SDSSr', 'CDS/P/HST/SDSSz']:
        url = (f"https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
               f"?hips={hips}&width={size}&height={size}"
               f"&fov={fov_deg}&projection=TAN&coordsys=icrs"
               f"&ra={ra}&dec={dec}&format=jpg")
        try:
            r = requests.get(url, timeout=15)
            if r.status_code == 200 and len(r.content) > 5000:
                return r.content
        except Exception:
            continue
    return None


def get_jwst_cutout(ra, dec, fov_deg=0.003, size=256):
    """Fetch JWST cutout via CDS HiPS2FITS, trying multiple band maps."""
    for hips in ['ESAVO/P/JWST/NIRCam_Imaging', 'CDS/P/JWST/F444W',
                 'CDS/P/JWST/F200W', 'CDS/P/JWST/F150W',
                 'CDS/P/JWST/F115W', 'CDS/P/JWST/EPO']:
        url = (f"https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
               f"?hips={hips}&width={size}&height={size}"
               f"&fov={fov_deg}&projection=TAN&coordsys=icrs"
               f"&ra={ra}&dec={dec}&format=jpg")
        try:
            r = requests.get(url, timeout=15)
            if r.status_code == 200 and len(r.content) > 5000:
                return r.content
        except Exception:
            continue
    return None


def get_euclid_cutout(ra, dec, radius_arcmin=0.1):
    """Fetch Euclid cutout. Try color (VIS+J+H), fall back to VIS grayscale."""
    coord = SkyCoord(ra, dec, unit='deg')
    rad_deg = radius_arcmin / 60.0

    query = (f"SELECT file_path, filter_name FROM q1.mosaic_product "
             f"WHERE INTERSECTS(CIRCLE('ICRS',{ra},{dec},{rad_deg}), fov)=1")
    results = Euclid.launch_job(query).get_results()
    if len(results) == 0:
        return None, None

    filters = {row['filter_name']: row['file_path'] for row in results}

    def fetch_array(fpath):
        tmp = os.path.join(os.getcwd(), '_euclid_cutout.fits')
        Euclid.get_cutout(file_path=fpath, coordinate=coord,
                         radius=radius_arcmin * u.arcmin, output_file=tmp)
        data = fits.getdata(tmp).astype(float)
        os.unlink(tmp)
        return data

    def norm(arr):
        p1, p99 = np.percentile(arr, [1, 99])
        return np.clip((arr - p1) / max(p99 - p1, 1e-10), 0, 1)

    def to_png(pil_img):
        buf = BytesIO()
        pil_img.save(buf, format='PNG')
        return buf.getvalue()

    # Try color: H=red, J=green, VIS=blue
    if all(k in filters for k in ['VIS', 'J', 'H']):
        try:
            arrays = {k: fetch_array(filters[k]) for k in ['VIS', 'J', 'H']}
            target = arrays['H'].shape
            for k in arrays:
                if arrays[k].shape != target:
                    arrays[k] = np.array(Image.fromarray(
                        arrays[k].astype(np.float32), mode='F'
                    ).resize((target[1], target[0]), Image.LANCZOS))
            rgb = np.stack([norm(arrays['H']), norm(arrays['J']),
                           norm(arrays['VIS'])], axis=-1)
            return to_png(Image.fromarray((rgb * 255).astype(np.uint8))), 'VIS+J+H'
        except Exception:
            pass

    for band in ['VIS', 'Y', 'J', 'H']:
        if band in filters:
            try:
                data = norm(fetch_array(filters[band]))
                return to_png(Image.fromarray((data * 255).astype(np.uint8))), band
            except Exception:
                continue

    return None, None


# --- Catalog queries -------------------------------------------------------

def ps1_gr_color(ra, dec, radius_arcsec=2.0):
    """Query PS1 DR2 for nearest source g-r color via MAST TAP."""
    radius_deg = radius_arcsec / 3600.0
    query = (f"SELECT TOP 1 gmeanpsfmag, rmeanpsfmag FROM dbo.meanobjectview "
             f"WHERE 1=CONTAINS(POINT('ICRS', ramean, decmean), "
             f"CIRCLE('ICRS', {ra}, {dec}, {radius_deg}))")
    url = 'https://mast.stsci.edu/vo-tap/api/v0.1/ps1dr2/sync'
    r = requests.get(url, params={'REQUEST': 'doQuery', 'LANG': 'ADQL',
                                  'FORMAT': 'csv', 'QUERY': query}, timeout=30)
    r.raise_for_status()
    lines = r.text.strip().split('\n')
    if len(lines) < 2:
        return None, None, None
    vals = lines[1].split(',')
    g = float(vals[0]) if vals[0] and float(vals[0]) != -999.0 else None
    r_mag = float(vals[1]) if vals[1] and float(vals[1]) != -999.0 else None
    gr = round(g - r_mag, 2) if g and r_mag else None
    return g, r_mag, gr


def skymapper_colors(ra, dec, radius_arcsec=2.0):
    """Query SkyMapper DR4 for nearest source g, r, v PSF mags via cone search."""
    radius_deg = radius_arcsec / 3600.0
    url = (f"https://skymapper.anu.edu.au/sm-cone/public/query"
           f"?RA={ra}&DEC={dec}&SR={radius_deg}&RESPONSEFORMAT=CSV&VERB=3")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    lines = r.text.strip().split('\n')
    if len(lines) < 2:
        return None, None, None, None, None
    header = lines[0].split(',')
    vals = lines[1].split(',')
    row = dict(zip(header, vals))
    def get_mag(key):
        v = row.get(key, '')
        return float(v) if v and v != 'NaN' else None
    g = get_mag('g_psf')
    r_mag = get_mag('r_psf')
    v = get_mag('v_psf')
    gr = round(g - r_mag, 2) if g and r_mag else None
    vg = round(v - g, 2) if v and g else None
    return g, r_mag, v, gr, vg


def hst_jwst_coverage(ra, dec, radius_arcsec=6.0):
    """Check for HST/JWST imaging at position. Returns dict of collections."""
    obs = Observations.query_criteria(
        coordinates=f"{ra} {dec}", radius=radius_arcsec / 3600.0,
        obs_collection=["HST", "JWST"], dataproduct_type="image")
    if len(obs) == 0:
        return {}
    result = {}
    for coll in set(obs["obs_collection"]):
        mask = obs["obs_collection"] == coll
        filters = sorted(set(obs["filters"][mask]))
        result[coll] = {"n_obs": int(mask.sum()), "filters": filters}
    return result


# --- Utilities -------------------------------------------------------------

def nJy_to_AB(flux_nJy):
    """Convert flux in nJy to AB magnitude."""
    if flux_nJy is None or flux_nJy <= 0:
        return None
    return -2.5 * np.log10(flux_nJy) + 31.4


def fmt_sigma(val, err):
    """Format value +/- error with sigma = |val/err|."""
    if val is None or err is None or err == 0:
        return '-'
    sigma = abs(val / err)
    return f'{val:.2f} +/- {err:.2f} ({sigma:.1f}s)'


def add_scale_bar(img_bytes, pixscale_arcsec, bar_arcsec=1.0, fmt='jpeg'):
    """Draw a horizontal 1 arcsec scale bar on the lower-right of a cutout image."""
    img = Image.open(BytesIO(img_bytes)).convert('RGB')
    draw = ImageDraw.Draw(img)
    bar_px = int(round(bar_arcsec / pixscale_arcsec))
    margin = 5
    y = img.height - margin
    x_right = img.width - margin
    x_left = x_right - bar_px
    draw.line([(x_left, y), (x_right, y)], fill='white', width=1)
    buf = BytesIO()
    img.save(buf, format=fmt.upper())
    return buf.getvalue()


# --- ANTARES search --------------------------------------------------------

STELLAR_CACHE = "stellar_ids.txt"


def load_stellar_ids(path=STELLAR_CACHE):
    """Load cached stellar locus IDs (one per line)."""
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def save_stellar_id(locus_id, path=STELLAR_CACHE):
    """Append one stellar locus ID to cache file."""
    with open(path, 'a') as f:
        f.write(locus_id + '\n')


def search_by_tag(tag, limit=10, skip_stellar=True):
    """Search ANTARES for loci with a given tag, in random order.
    If skip_stellar=True, filters out Gaia >3sigma PM/parallax on the fly,
    caching stellar IDs so re-runs skip them instantly.
    """
    query = {
        "query": {
            "function_score": {
                "query": {"bool": {"filter": [{"term": {"tags": tag}}]}},
                "random_score": {}
            }
        }
    }
    stellar_ids = load_stellar_ids() if skip_stellar else set()
    results = []
    n_cached = 0
    n_checked = 0
    n_total = 0
    max_total = limit * 10  # stop after checking this many loci total
    for locus in search(query):
        n_total += 1
        if n_total > max_total:
            print(f"  (reached {max_total} loci limit, stopping)")
            break
        if locus.locus_id in stellar_ids:
            n_cached += 1
            if n_cached % 100 == 0:
                print(f"  ... {n_cached} cached stellar skipped so far")
            continue
        if skip_stellar and is_stellar(locus):
            n_checked += 1
            save_stellar_id(locus.locus_id)
            stellar_ids.add(locus.locus_id)
            if n_checked % 100 == 0:
                print(f"  ... {n_checked} new stellar found so far")
            continue
        results.append(locus)
        if len(results) >= limit:
            break
    print(f"  checked {n_total} loci: {len(results)} kept, "
          f"{n_cached} cached stellar, {n_checked} new stellar")
    return results


def is_stellar(locus, sigma=3.0):
    """Check if locus has Gaia PM or parallax > sigma-sigma (stellar contaminant).
    Returns True if any of pmra, pmdec, or parallax exceeds threshold.
    """
    cats = list(locus.catalogs or [])
    if "gaia_dr3_gaia_source" not in cats:
        return False  # no Gaia match -> not confirmed stellar
    gaia_rows = locus.catalog_objects.get("gaia_dr3_gaia_source", [])
    for row in gaia_rows:
        for key in ["pmra", "pmdec", "parallax"]:
            val = row.get(key)
            err = row.get(f"{key}_error")
            if val is not None and err is not None and err > 0:
                if abs(val / err) > sigma:
                    return True
    return False


# --- CSV classification I/O ------------------------------------------------

CSV_COLUMNS = ["locus_id", "ra", "dec", "grade", "reason", "note",
               "ps1_g_r", "gaia_bp_rp", "sci_mag", "tmpl_mag", "timestamp"]


def save_to_csv(csv_path, locus_id, ra, dec, grade, reason="", note="",
                ps1_g_r="", gaia_bp_rp="", sci_mag="", tmpl_mag=""):
    """Append one classification row to CSV. Creates file with header if needed."""
    write_header = not os.path.exists(csv_path)
    with open(csv_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(CSV_COLUMNS)
        writer.writerow([locus_id, ra, dec, grade, reason, note,
                         ps1_g_r, gaia_bp_rp, sci_mag, tmpl_mag,
                         datetime.now().strftime("%Y-%m-%d %H:%M:%S")])


def load_classified_ids(csv_path):
    """Return set of locus_ids already classified in the CSV."""
    if not os.path.exists(csv_path):
        return set()
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        return {row["locus_id"] for row in reader}
