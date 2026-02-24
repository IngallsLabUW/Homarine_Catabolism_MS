# Homarine_Catabolism
Data analysis for homarine catabolism manuscript.  Preprint at https://www.researchsquare.com/article/rs-7359689/v1.

## Data Retrieval

Much of the data needed are not hosted on this github due to sizes of data.  To retrieve raw data for data analysis, you will need to populate data folders using scripts in the data folder. See **data/README.md** for more information.

## Repository Organiziation

### **data/raw**
Raw data files (unedited as much as possible).  

### **data/intermediate**
Intermediate data products. 

### **R_functions**
R functions used across multiple scripts.  Separated into logical separate .R scripts and documented using Roxygen documentation.  


## Dependency Management

### R

This project uses `renv` for package management.  After cloning the github repository, open the R project and run `renv::restore()` to make sure your packages match. To learn more about how renv works, [see this resource](https://rstudio.github.io/renv/articles/renv.html).

### Python

This project uses pip paired with venv to manage dependencies. Note that requirements_dev.txt should be used for development dependencies, and requirements.txt should be used for production/binder dependencies (added manually and with discretion).

#### To install the dependencies:

1. Clone the github repository
2. create a virtual environment:
    `python -m venv venv`
3. Activate the virtual environment:
    `source venv/bin/activate`
4. Install the necessary packages:
    `pip install -r requirements.txt`
    **Note** to update your package installations:
        `pip install -U -r requirements.txt`

#### To add new packages:

1. Activate the virtual environment:
    `source venv/bin/activate`
2. Install any new packages:
    `pip install <package>`
3. Capture the new requirements:
    `pip freeze > requirements_dev.txt`
4. Push changes to github




