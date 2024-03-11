# Readme for raw data and associated data retrieval scripts

## Metabolomics data overview

### Obi1 metabolomics data
Data for the Obi1 metabolomics experiment are populated by sourcing the *data/retrieve_metabolomics_data.R* script from the project level directory.  

Note that this script is a bit slow to fully download all the data. If you already have all the data, this can be skipped. If the script is run and you already have all the data, no worries! It will skip over data you already have downloaded.

This will create and populate the **data/raw/metabolmics/obi1** folder and its subfolders.