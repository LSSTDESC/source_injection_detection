"""
Fix point source RA/DEC metadata in H5 files using CSV as ground truth.

Problem: In some H5 files, metadata.attrs['point_source_light_{band}_ra_image_{ind2}']
and the corresponding dec values are incorrect due to ID misalignment between CSV and H5.

Specifically, file 1 has a missing ID (col0=24) in the CSV, so H5 indices >= 24 are
offset by 1 relative to CSV col0. Files 2-7 have no missing IDs.

This script reads the correct values from the CSV files and updates the H5 metadata
in-place. Run this BEFORE merge_h5.py and merge_clean_h5.py.

Usage:
    python fix_h5_metadata.py
"""

import numpy as np
from astropy.table import Table

# h5py may not be available in all environments
import h5py

BANDS = ['u', 'g', 'r', 'i', 'z', 'y']
N_FILES = 7
FOLDER = "lens_finding_postage_stamps_dataset"


def get_csv_id_mapping(setid):
    """Return a mapping from H5 index to CSV row for a given setid.

    For file 1: col0=24 is missing, so H5 ind >= 24 maps to CSV col0 = ind + 1.
    For files 2-7: direct mapping (H5 ind = CSV col0).
    """
    csv = Table.read(f"{FOLDER}/3000sqdeg_lsst_1y_sample_{setid}.csv")
    col0 = np.array(csv['col0'])

    # Find missing IDs
    all_ids = set(range(max(col0) + 1))
    missing = sorted(all_ids - set(col0))

    # Build lookup: col0 value -> row index in csv table
    col0_to_row = {c: i for i, c in enumerate(col0)}

    return csv, col0_to_row, missing


def fix_file(setid):
    """Fix point source metadata in one H5 file using its corresponding CSV."""
    csv, col0_to_row, missing = get_csv_id_mapping(setid)

    h5_filename = f"{FOLDER}/3000sqdeg_lsst_1y_sample_{setid}.h5"
    print(f"\n{'='*60}")
    print(f"Processing setid={setid}: {h5_filename}")
    print(f"  CSV rows: {len(csv)}, missing IDs: {missing}")

    n_updated = 0
    n_mismatch = 0

    with h5py.File(h5_filename, 'r+') as hf:
        n_systems = len(hf.keys())
        print(f"  H5 systems: {n_systems}")

        for ind in range(n_systems):
            key = f"lsst_lens_{ind}"
            if key not in hf:
                continue

            # Map H5 index to CSV col0
            csv_col0 = ind
            for m in missing:
                if ind >= m:
                    csv_col0 += 1

            if csv_col0 not in col0_to_row:
                print(f"  WARNING: H5 {key} maps to col0={csv_col0} which is not in CSV, skipping")
                continue

            row_idx = col0_to_row[csv_col0]
            row = csv[row_idx]

            metadata = hf[key]['metadata']
            light_curves = hf[key]['light_curves']
            n_images = len(light_curves)

            for band in BANDS:
                for ind2 in range(n_images):
                    ra_key = f'point_source_light_{band}_ra_image_{ind2}'
                    dec_key = f'point_source_light_{band}_dec_image_{ind2}'

                    csv_ra = row[ra_key]
                    csv_dec = row[dec_key]

                    h5_ra = metadata.attrs[ra_key]
                    h5_dec = metadata.attrs[dec_key]

                    if not np.isclose(h5_ra, csv_ra, atol=1e-10) or \
                       not np.isclose(h5_dec, csv_dec, atol=1e-10):
                        n_mismatch += 1
                        if n_mismatch <= 5:
                            print(f"  Mismatch: {key} {ra_key}: H5={h5_ra:.6f} CSV={csv_ra:.6f}")

                    metadata.attrs[ra_key] = csv_ra
                    metadata.attrs[dec_key] = csv_dec
                    n_updated += 1

    print(f"  Updated {n_updated} attrs, found {n_mismatch} mismatches")


def main():
    for setid in range(1, N_FILES + 1):
        fix_file(setid)
    print(f"\nDone. All {N_FILES} H5 files updated.")


if __name__ == "__main__":
    main()
