# Homarine_Catabolism
Data analysis for homarine catabolism manuscript

## Dependency Management
This project uses `renv` for package management.  After cloning the github repository, open the R project and run `renv::restore()` to make sure your packages match. To learn more about how renv works, [see this resource](https://rstudio.github.io/renv/articles/renv.html).

## Data Retrieval
**Work in progress not yet functional**

Much of the data needed are not hosted on this github due to sizes of data.  To retrieve raw data for data analysis, you will need to populate data folders using scripts in the data folder. See **data/README.md** for more information.

## Repository Organiziation

### **data/raw**
Raw data files (unedited as much as possible).  

### **data/intermediate**
Intermediate data products. 

### **R_functions**
R functions used across multiple scripts.  Separated into logical separate .R scripts and documented using Roxygen documentation.  




