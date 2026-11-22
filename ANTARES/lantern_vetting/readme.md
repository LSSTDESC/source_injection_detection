`lantern_viewer.ipynb` and `lantern_viewer_v1.ipynb` check Legacy Survey images. These will be replaced by checking DP2 deep coadd images.

`dp2_deep_coadd*.ipynb`: Run on RSP, see DP2 deep coadd --- `full` can see the whole patch, `cutout` just gets a stamp without meta data, `cutout-exposure` includes meta data (PSF model etc).

`lantern_candidate_rgb_vetting.ipynb`: After second pass and cleaning, looking at the DP2 coadd images.

---

`vetting.ipynb` checks a single locus, while `vetting_2.ipynb` goes through a few loci and takes notes (more details in `vetting_lib.py`).

## Updated (uploaded Nov 2026)

`lantern_candidate_rgb_vetting_new.ipynb`: take the catalog of targets after 2nd pass and Gaia pm/parallax cuts; early vetting (saves PDF rather than HTML, see below); extendedness and PSF–aperture flux cuts; known-lens comparisons

`lantern_candidate_cutout_viewer.ipynb`: visual inspection using DP2, Legacy and ANTARES images; generates PDFs
