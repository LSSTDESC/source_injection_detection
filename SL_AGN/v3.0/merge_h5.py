import h5py
import re
from pathlib import Path

INPUT_FILES = [f"3000sqdeg_lsst_1y_sample_{i}.h5" for i in range(1, 8)]
OUTPUT_FILE = "3000sqdeg_lsst_1y_sample_merged.h5"
INPUT_DIR = Path(".")
KEY_PREFIX = "lsst_lens"
BANDS = ['u', 'g', 'r', 'i', 'z', 'y']

KEY_PATTERN = re.compile(rf"{re.escape(KEY_PREFIX)}_(\d+)$")

def key_sort_fn(name):
    m = KEY_PATTERN.match(name)
    # get the integer part
    return int(m.group(1)) if m else float("inf")

def is_valid_system(hf, key):
    """Check if a lens system has complete data across all 6 bands."""
    system = hf[key]

    # Check required top-level groups exist
    for group in ['static_image', 'light_curves', 'postage_stamps', 'observation_dates']:
        if group not in system:
            return False

    # Check 'image_0' exists under light_curves
    if 'image_0' not in system['light_curves']:
        return False

    for band in BANDS:
        # Check static image exists
        if band not in system['static_image']:
            return False
        if 'lens_plus_lensed_agn_host' not in system['static_image'][band]:
            return False

        # Check light curve exists
        if band not in system['light_curves']['image_0']:
            return False

        # Check postage stamps exist
        if band not in system['postage_stamps']:
            return False
        if 'all_observations' not in system['postage_stamps'][band]:
            return False
        n_stamps = system['postage_stamps'][band]['all_observations'].shape[0]
        if n_stamps < 1:
            return False

        # Check observation dates exist and are consistent
        if band not in system['observation_dates']:
            return False
        n_dates = len(system['observation_dates'][band])
        if n_dates < 1 or n_stamps != n_dates:
            return False

        # Epoch limit (consistent with stamp.py)
        if n_stamps > 30:
            return False

    return True

def main():
    global_idx = 0
    total_scanned = 0
    total_skipped = 0
    with h5py.File(OUTPUT_FILE, "w") as fout:
        for fname in INPUT_FILES:
            fpath = INPUT_DIR / fname
            if not fpath.exists():
                raise FileNotFoundError(f"Missing input: {fpath}")

            with h5py.File(fpath, "r") as fin:
                # sort by the integer part of the key name
                matching_keys = sorted(
                    [k for k in fin.keys() if KEY_PATTERN.match(k)],
                    key=key_sort_fn,
                )
                file_total = len(matching_keys)
                file_skipped = 0
                print(f"\n{fname}: {file_total} groups")

                for old_key in matching_keys:
                    total_scanned += 1
                    if not is_valid_system(fin, old_key):
                        file_skipped += 1
                        total_skipped += 1
                        continue
                    new_key = f"{KEY_PREFIX}_{global_idx}"
                    fin.copy(fin[old_key], fout, name=new_key)
                    global_idx += 1

                print(f"  kept: {file_total - file_skipped}, skipped: {file_skipped}")

    print(f"\nDone. Scanned {total_scanned}, skipped {total_skipped}, wrote {global_idx} groups to {OUTPUT_FILE}.")

if __name__ == "__main__":
    main()
