Here we store code for downloading tagged alerts from ANTARES. 

The v1 API is based on the search function [here](https://nsf-noirlab.gitlab.io/csdc/antares/client/_modules/antares_client/search.html). 

The v2 API is listed [here](https://api.antares.noirlab.edu/docs).

## What has been updated

`query_v1.ipynb`: Basic query (1k loci).

`query_v2.ipynb`: Run query in parallel.

`query_v3.ipynb`: Add a few Gaia quantities;

`query_v3-1.ipynb`: some functions (e.g. retry) for network;

`query_v3-2.ipynb`: Add a few Gaia quantities (10k loci).

`query_v4.ipnb` (in preparation): new API?

---

`check_dates.ipynb`: make histograms for the alert dates (based on the locus query result).

`skip_alerts_before_May_27.ipynb`: skip the alerts before the `lantern` filter was implemented; clean empty loci; remove duplicated loci and alerts.

`check_info_visualization.ipynb`: make sky plots for the cleaned loci/alerts.

`lantern_second_pass_local_candidate.ipynb`: run 2nd pass on the cleaned loci/alerts to get candidates.

`lantern_second_pass_local_candidate_multi_loci.ipynb`: run 2nd pass on the cleaned loci/alerts to get candidates; if some loci are close, treat as one target.

`lantern_candidate_analysis_after_2nd_pass.ipynb`: make plots for the candidates after the 2nd pass; select new candidates by Gaia info, mark milliquas.



---

`API-v2-endpoint.ipynb`: initial test using the ANTARES' API v2.
