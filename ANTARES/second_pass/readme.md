`lantern_second_pass_local_candidate.ipynb`: run 2nd pass on the cleaned loci/alerts to get candidates.

`lantern_second_pass_local_candidate_multi_loci.ipynb`: run 2nd pass on the cleaned loci/alerts to get candidates; if some loci are close, treat as one target.

`lantern_candidate_analysis_after_2nd_pass.ipynb`: make plots for the candidates after the 2nd pass; select new candidates by Gaia info, mark milliquas.


## Updated

`lantern_second_pass_multi_loci.ipynb`: run second pass on the downloaded loci (from query_v4; group loci within 3 arcseconds, load second-pass pickled model, calculate second-pass score, and select candidates).

`lantern_candidate_analysis_multi_loci.ipynb`: analyze second-pass candidates, calculate photometric and Gaia properties, remove likely stars using significant Gaia proper motion or parallax (>3σ).

`lantern_candidate_analysis_multi_loci_extra_plots.ipynb`: supplementary notebook for `lantern_candidate_analysis_multi_loci.ipynb` showing extra diagnostic plots.

`lantern_variable_star_crossmatch`: version a and b; take input `targets.csv` with columns "name,ra,dec" and check variable stars in Gaia, AAVSO, OGLE (version b omits Gaia due to speed).
