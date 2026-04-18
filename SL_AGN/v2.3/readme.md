Following v2.2, this version implements item 3 of the v2.2 plan:

- More epochs
    - Pick a random visit per run (from the valid range in `BAND_EXPOSURE_TOTAL` in `lib/tools_template.py`)
    - Extract all available observation epochs per stamp system at step0
    - Randomly select one epoch per stamp at step1 injection
