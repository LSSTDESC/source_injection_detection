# source_injection_detection

This repo is for storing code for source injection and image reprocessing, primarily on Rubin images and using the LSST Science Pipelines.  

Example DESC projects using this repo include:
* "The First Cut: Developing the Initial Catalog Query for LSST Lensed AGNs".

The idea is to collect notebooks and scripts in one place so they can be picked up, adapted and re-used in new science cases.

Organization is roughly one folder per science case.
* `SL_AGN` This folder stores the code for the source injection of strongly lensed AGN in LSST. Go to the folder for one version of the code and test it (`bash run.sh`). 

* `ANTARES` This folder stores the code for studying ANTARES alerts and ZTF data for lens systems and candidates.

* `DP1` This folder stores the code for some basic analysis of the DP1 data.

Note: the code for lensed AGN search -- making broker filters, analyzing simulated DIA sources, analyzing clusters, is stored in `lantern` [here](https://github.com/drphilmarshall/lantern/tree/main).
