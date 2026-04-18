import h5py
import re
from pathlib import Path

INPUT_FILES = [f"3000sqdeg_lsst_1y_sample_{i}.h5" for i in range(1, 8)]
OUTPUT_FILE = "3000sqdeg_lsst_1y_sample_merged.h5"
INPUT_DIR = Path(".")
KEY_PREFIX = "lsst_lens"

KEY_PATTERN = re.compile(rf"{re.escape(KEY_PREFIX)}_(\d+)$")

def key_sort_fn(name):
    m = KEY_PATTERN.match(name)
    # get the integer part
    return int(m.group(1)) if m else float("inf")

def main():
    global_idx = 0
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
                print(f"\n{fname}: {len(matching_keys)} groups ")
                print("sorted keys: ", matching_keys)

                for old_key in matching_keys:
                    new_key = f"{KEY_PREFIX}_{global_idx}"
                    fin.copy(fin[old_key], fout, name=new_key)
                    global_idx += 1

    print(f"\nDone. Wrote {global_idx} groups to {OUTPUT_FILE}.")

if __name__ == "__main__":
    main()
