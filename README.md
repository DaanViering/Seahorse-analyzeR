Seahorse-analyzeR is an R-based programme developed to analyze Seahorse XF96e data.
More information can be found in the scientific publication that will follow soon.
Please contact me if you have any questions or requests for additional features.

The script_main2.R can be run to analyze Seahorse XF96e data. 

The current version only supports analysis of the maximal respiratory capacity.
The scripts in add_columns.R, configure.R and load_functions2.R are called upon by script_main2.R .
bind_all_rate_data.R can be used to combine rate data files from multiple Seahorse experiments into one file for direct analysis. To use it, place all rate data files in a folder with the name 'rate_data', directly in the parent folder 'Seahorse-analyzeR'.
normalization_maximal_resp.xlsx can be used to subsequently normalize the data obtained with Seahorse-analyzeR to all control cell lines. Future versions of Seahorse-analyzeR will incorporate this in the main R-script.

The script was written in R, version 3.6.2 for Mac, using Rstudio version 1.1.456. The R-script depends on packages data.table (version 1.12.8) and tidyverse (1.3.0).
