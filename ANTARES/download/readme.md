Here we store code for downloading tagged alerts from ANTARES. 

The v1 API is based on the search function [here](https://nsf-noirlab.gitlab.io/csdc/antares/client/_modules/antares_client/search.html). 

The v2 API is listed [here](https://api.antares.noirlab.edu/docs).

## What has been updated

`query_v1.ipynb`: Basic query (1k loci).

`query_v2.ipynb`: Run query in parallel.

`query_v3.ipynb`: Add a few Gaia quantities;

`query_v3-1.ipynb`: some functions (e.g. retry) for network;

`query_v3-2.ipynb`: Add a few Gaia quantities (10k loci).

## Updated (uploaded Nov 2026)

`query_v4.ipnb`: new API (v2) on ANTARES for loci/alerts. There are two v4 notebooks currently (download 10k and (additional) 100k loci). Elasticsearch has a limit of 10k per query (result ordering unchanged by default), and the 2nd one use a random score (seed) to obtain different samples. The 2nd notebook also gets Gaia info via v1 API (so it's "hybrid"; this can be updated). In the future we can try downloading loci by RA; there can still be a 10k limit, which can be solved by splitting the angle range by half each time.

`antares_loci_alerts_visualization.ipynb`: visualizing the downloaded loci/alerts.

## Others

`check_dates.ipynb`: make histograms for the alert dates (based on the locus query result).

`skip_alerts_before_May_27.ipynb`: skip the alerts before the `lantern` filter was implemented; clean empty loci; remove duplicated loci and alerts.

`check_info_visualization.ipynb`: make sky plots for the cleaned loci/alerts.



---

`API-v2-endpoint.ipynb`: initial test using the ANTARES' API v2.
