Note: in order to generate datasets, you **need to have GPU access**.

You can generate a dataset by calling:
``python generate_slsim_lsst_dataset.py --sky-area 10000 --n-iter 1 --baseline 10 --output-file data/test.h5``

This will sample over 10000 sq deg once and create time-series postage stamp images of lensed AGN over a 10 year baseline and save it to data/test.h5. Make sure the `data` directory exists.

You can also add the following flags to you call to customize the dataset selection function:
1. `--image-separation` (this is the lower bound of the maximum separation between two lensed images; default = 0.5 arcsec)
2. `--magnitude-limit` (this is the lensed magnitude of the second brightest image; default = 24 mag)
3. `--start-index` (default: 0)
4. `--end-index` (default: number of lenses in your dataset)
5. `--log-file` (place to log the processing and how many lenses were successfully stored, etc; default: lens_generation_log.json)

If you're running this on a cluster, you can also just submit this as a job using sbatch. It roughly takes ~1.5 hours to generate a 3000 sq deg 1-year baseline dataset.

Here's the `sbatch` call I use with `generate_slsim_auto.sh`: 
> sbatch generate_slsim_auto.sh 1000 3 1 3000sqdeg_lsst_1y_sample_5.h5

Note that you should edit the bash script with the prefix for your dataset target location (it is currently my scratch directory: `/pscratch/sd/v/vpadma/lens_finding/` )