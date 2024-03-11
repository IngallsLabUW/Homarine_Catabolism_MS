# Homarine_Catabolism
Data analysis for homarine catabolism manuscript

## Dependency Management
This project uses `renv` for package management.  After cloning the github repository, open the R project and run `renv::restore()` to make sure your packages match. To learn more about how renv works, [see this resource](https://rstudio.github.io/renv/articles/renv.html).

## Data Retrieval
**Work in progress not yet functional**
Much of the data needed are not hosted on this github due to sizes of data.  To retrieve data and populate data folders, run `source("data/retrieve_data.R")`. *Note that this script is a bit slow to fully download all the data*. If you already have all the data, this can be skipped.  If the script is run and you already have all the data, no worries!  It will skip over data you already have downloaded.

## Repository Organiziation

### **data/raw**
Raw data files (unedited as much as possible).  
See **data/raw/README.md** for more information.

### **data/intermediate**
Intermediate data products. 
See **data/intermediate/README.md** for more information.

### **R_functions**
R functions used across multiple scripts.  Separated into logical separate .R scripts and documented using Roxygen documentation.  




