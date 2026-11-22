`lantern_viewer.ipynb` and `lantern_viewer_v1.ipynb` check Legacy Survey images. These will be replaced by checking DP2 deep coadd images.

`dp2_deep_coadd*.ipynb`: Run on RSP, see DP2 deep coadd --- `full` can see the whole patch, `cutout` just gets a stamp without meta data, `cutout-exposure` includes meta data (PSF model etc).

`lantern_candidate_rgb_vetting.ipynb`: After second pass and cleaning, looking at the DP2 coadd images.

---

`vetting.ipynb` checks a single locus, while `vetting_2.ipynb` goes through a few loci and takes notes (more details in `vetting_lib.py`).

## Updated (uploaded Nov 2026)

`lantern_candidate_rgb_vetting_new.ipynb`: take the catalog of targets after 2nd pass and Gaia pm/parallax cuts; early vetting (PDF rather than HTML below); extendedness and PSF–aperture flux cuts (looks for extended ones rather than point sources; can be adjusted for point sources or extended ones in DIA); known-lens comparisons.

`lantern_candidate_cutout_viewer.ipynb`: visual inspection using DP2, Legacy and ANTARES images; generates PDFs (A similar one is `lantern_candidate_cutout_viewer_red.ipynb` which checks red objects); takes catalogs from `lantern_candidate_rgb_vetting_new.ipynb`.

`lantern_individual_target_html_viewer`: (notebook or python script) detailed individual-target HTML, including coadds, alerts, PSFs and candidate metadata.

`lantern_dp2_irg_candidate_grid.ipynb`: DP2 cutouts; `lantern_dp2_predicted_irg_grid.ipynb`: SharPy deconvolution; `predicted_irg_grid_auto_per_image_band`: automatic linear min/max instead of fixed min/max.

`lens_cutouts_space_telescope.ipynb`: cutout for space telescope images.

`dp2_position_angle.ipynb`: measure the position angle of a point source pair.
