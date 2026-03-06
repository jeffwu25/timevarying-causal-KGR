# Timevarying-Causal-KGR

This repository contains the code and data used in the paper: 

“Causal effects of time-varying wildfire smoke exposure on respiratory-related mortality in California.”

Authors: Authors: Jeffrey Wu, Gareth W Peters, Alex Franks (UCSB)

In this project, we estimate the causal effect of monthly wildfire specific PM2.5 accumulation (smoke PM2.5) on respiratory-related mortality in California from 2014-2019. While most existing studies evaluate daily or annual smoke exposure, this work focuses on cumulative smoke PM2.5 exposure over several months, which may better capture delayed health effects from prolonged wildfire events. 

To estimate causal effects, we account for time-varying confounding and treatment-confounder feedback using marginal structural models (MSMs) with stabilized inverse probability weighting (IPW). We estimate a nonlinear, causal average dose-response function (ADRF) relating cumulative smoke PM2.5 to respiratory-related mortality. To investigate heterogeneity in health effects across socioeconomic contexts using social deprivation index (SDI). We model our exposure and outcome processes using kernel graph regression (KGR), a Bayesian spatiotemporal modeling framework that captures complex spatiotemporal dependence structures in a parsimonious and interpretable fashion. 

This framework provides a basis for evaluating spatiotemporal causal inference problems with a continuous, time-varying treatment/exposure. The high-level pipeline is as follows: 
Data -> Estimate spatial and temporal dependence -> Estimate stabilized IPWs with exposure model -> Estimate dose-response function with outcome model -> Estimate and visualize causal effects 


## Repository Structure: 
```
Timevarying-Causal-KGR/
│
├── Data/
│
├── Exploratory Analysis/
│
├── Main Analysis/
│
└── Sample Output/
```

## Data Folder

The Data folder contains all the datasets and geographic files needed to run the analysis code in P2-Main-Analysis-cum2 (and cum6). These datasets were downloaded from various publicly available sources. 

- CA_census_pops1019.xlsx: contains census population data for counties of California for each year 2010-2019
- CA_shapefiles.zip: contains shape files needed to perform spatial analysis and plot maps of California. Shapefiles are a common format for vector-based geographic information system (GIS) data. They can be opened and used in any GIS software and in R or Python. A shapefile consists of multiple file types beyond the .shp (specifically, .cpg, .dbf, .prj, .sbn, and .sbx). The user only interacts directly with the .shp file but the other files need to be in the same directory.
- Population_Categories.xlsx: contains annual estimates of the resident population for California state for selected age groups and by sex for the years 2020-2022
- SoA.data.1019.xlsx: contains county-level social deprivation scores and associated subindices from a report by the Society of Actuaries
- Cal-ViDa_Death1423.xlsx: contains respiratory-related mortality counts queried from the California Department of Public Health’s Cal-ViDa query tool
- final_EPA_data.csv: contains air pollutant measurements queried and aggregated from the EPA Air Quality System AQI
  - Lead (micrograms/m^3)
  - CO (parts per million)
  - SO2 (parts per billion)
  - NO2 (parts per billion)
  - O3 (parts per million)
  - PM10 (micrograms/m^3)
  - PM2.5 (micrograms/m^3)
- cimis-monthly-weather.csv: contains meteorological measurements queried and aggregated from the California Irrigation Management Information System (CIMIS) database
  - Precipitation (inches)
  - Temperature (Fahrenheit)
  - Wind (mph)
  - Humidity/dew point (Fahrenheit)
- flu-ili-byregion-fluseason.csv: contains measurement of flu and influenza-like-illness prevalence by region during flu season in California 
- final_inla_df.csv: contains all the monthly, county data (exposure, outcome, confounders) used for the main analysis
- wildfirepm25_monthly_total: contains monthly smoke PM2.5 levels for each county in California from 2014-2019


## Exploratory Analysis 

The Exploratory Analysis folder contains scripts that download, clean, and explore the data used in the main analysis

- P2-EDA.qmd: contains the code to
  - Download the Society of Actuaries’ SDI data -> explore trends by county and over time
  - Creating clusters of spatially contiguous counties based on SDI
  - Download the EPA’s air pollutant data and ensure stationarity for each time series (detailed procedure in Download-EPA.qmd)
  - Download the California Department of Public Health’s Cal-ViDa mortality data and impute censored values with an EM algorithm before aggregating to the monthly, county level -> explore spatial and temporal patterns
  - Download Childs et. al.’s smoke PM2.5 data and aggregate it to the monthly (from daily), county level -> compare with recent notable wildfire events to check the accuracy of their predictions
  - Download CIMIS’ station data and construct monthly summaries for each weather variable for each county
  - Download the California Department of Public Health’s ILI-prevalence data and construct monthly county level prevalence summaries from the regional surveys
  - Join these datasets together into a single dataset (inla_df_final) for the main analysis
  - Group counties into quintile groups based on SDI
  - Plot time series for all variables to ensure stationarity and minimal measurement error
  - Analyze the residual component of each time series and investigate potential exposure-confounder feedback with cross-lag autocorrelations, Spearman correlation, and Granger causality

- Download-EPA.qmd: contains code to download station data from the EPA AQS and then aggregate it into monthly, county level summaries for each air pollutant. This involved:
  - Querying all stations in California that measure at least one of the air pollutants

Then, for each air pollutant and each county, we

  - Identify a subset of stations that are close enough to provide a reasonable spatial coverage
  - Filter out any stations with too many missing values or outliers
  - Combine the data from the different stations into a monthly median (except for ambient PM2.5 because we wanted to have a monthly total)

This produces a dataset (final_EPA_data.csv) measuring each air pollutant’s monthly median concentration in each county between 2014-2019. We directly import this data into P2-EDA.qmd. We perform a similar procedure to obtain county-level monthly summaries from the CIMIS station data in P2-EDA.qmd.


## Main Analysis

The Main Analysis folder contains scripts that carry out the proposed causal inference analysis and visualize the results 

- P2-Main-Analysis-cum2.qmd: evaluates the causal effect of a two month accumulation of smoke PM2.5 on respiratory-related mortality in California from 2014-2019. Important steps include:
  - Creating lagged and cumulative time series for certain variables, most notably smoke PM2.5
  - Estimating a spatial (with a graph filter) and temporal dependence structure (with a kernel)
  - Estimating a Log Normal exposure model with INLA to estimate stabilized IPWs
  - Performing posterior predictive checks to validate the choice of a Log Normal model over a Gamma model
  - Estimating a Negative Binomial outcome model with INLA which is used to make counterfactual predictions and characterize the causal ADRF
  - Performing posterior predictive checks to validate the choice of a Negative Binomial model over a Poisson model
  - Evaluate out of sample prediction performance for exposure and outcome models
  - Use predicted counterfactual outcomes from Estimate-DRFs.R output to construct causal ADRF and estimate causal effects (population and conditional average treatment effects and risk ratios)
  - Present results with various visualizations and tables

- P2-Main-Analysis-cum6.qmd: evaluates the causal effect of a six month accumulation of smoke PM2.5 on respiratory-related mortality in California from 2014-2019
  - Same code/steps as cum2 analysis

- Estimate-DRFs.R: predicts counterfactual outcomes for all county-months for different levels of cumulative smoke PM2.5 and saves the results in a list (run as a background job). Important steps include:
  - Load workspace from main analysis (data, stabilized IPWs, etc)
  - Define a grid of exposure values based on quantiles of cumulative smoke PM2.5
  - For all county-month observations, make counterfactual predictions for a given level of exposure (how many respiratory-related deaths would we expect for a given level of exposure, for each county and month)
  - Store predictions for each dose grid level in a list and export
  - This takes about 8 hours to run on 384 cores because we make predictions for 58 counties (which have to be done in 4 batches) for 10 dose grid levels

- Estimate-DRFs-wCI.R: generates bootstrapped datasets with wild bootstrap and then estimates DRFs with each of the bootstrap samples to obtain uncertainty quantification. Important steps include:
  - Load workspace from main analysis (data, stabilized IPWs, etc)
  - Define a function to estimate DRFs (same one from Estimate-DRFs.R)
  - Define a function to perform wild bootstrap
  - Define a parallel computing structure to estimate 25 DRFs at a time
  - Define a grid of exposure values based on quantiles of cumulative smoke PM2.5 (same one from Estimate-DRFs.R)
  - Store estimated DRFs for each bootstrapped dataset in a list and export
  - This takes about 70 hours to run on 384 because we are running the Estimate-DRFs.R script on 25 different bootstrapped datasets in parallel. Running too many parallel jobs at once puts a significant strain on computational resources


## Sample Output

The Sample Output folder contains some of the significant results from our causal analysis in addition to a knitted HTML of our main analysis QMD. 

•	P2-Main-Analysis-cum2 HTML: knitted HTML of main analysis
•	CCF plots to demonstrate treatment confounder feedback: cross correlation function plots between ILI prevalence and the exposure/outcome  illustrates ILI prevalence is a potential confounder and there is potentially exposure-confounder feedback to account for 
•	ADRF plot with CI: causal ADRF estimated for cum2 smoke PM2.5
•	Relative risk ratio curve with CI: relative risk ratios calculated from the causal ADRF between exposure = 0 and every other value on the dose grid 
•	Relative risk ratio curve by SDI quintile group with CI: relative risk ratios calculated from the SDI quintile group specific causal ADRFs between exposure = 0 and every other value on the dose grid 


## Software Requirements: 

This analysis was run on a virtual computing server called Roble with 384 cores - running the analysis locally will be much slower. This analysis was conducted in R (v4.4.3) using the following packages: 

### Data manipulation and utilities
• dplyr  
• tidyverse  
• stringr  
• lubridate  
• zoo  
• reshape2  
• jsonlite  
• readxl  
• utils  
• pryr  
• tictoc  

### Visualization
• ggplot2  
• ggh4x  
• ggbeeswarm  
• patchwork  
• corrplot  
• superheat  
• viridis  
• RColorBrewer  
• scales  
• gt  
• webshot2  

### Spatial analysis
• sp  
• spdep  
• geosphere  
• urbnmapr  

### Bayesian and statistical modeling
• INLA  
• brinla  
• invgamma  
• HMMpa  
• huge  
• pscl  
• splines  
• shapes  
• fmesher  

### Air quality data
• RAQSAPI  
• con2aqi  

### Parallel and high-performance computing
• parallel  
• foreach  
• doParallel  
• future  
• future.apply  
• RhpcBLASctl  

### Development tools
• devtools







