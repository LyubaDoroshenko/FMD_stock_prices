# Functional Motif Discovery in Stock Market Prices
This repository contains the R code and data files used in the article [Functional motif discovery in stock market prices](https://dx.doi.org/10.2139/ssrn.4642040) by Marzia A. Cremona, Lyubov Doroshenko and  Federico Severino.

The folders `Simulation_study` and `Empirical_analysis` enclose codes of the simulation study and empirical analysis, respectively. 

## Simulation study

The folder `Simulation study` encloses the following files and folders:

- `simulated_length_300.R` - R code for simulating 20 curves of length 300 following the Mixed Autoregressive MAR(1,1) process described in the Section "Simulation study". The resulting simulated data are in the file `data_smoothed.RData`.

- `simulated_300_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated in the file `simulated_length_300.R`, using partially random initialisation and alpha=0.9 (see Section "Simulation study"). The directory `sim300_fmd` contains the candidate motifs resulting from ProbKMA-FMD. This allows to reduce the computational cost by directly uploading the candidate motifs without re-running all steps of ProbKMA-FMD. The code performs a new run of ProbKMA-FMD in the case the result directory `sim300_fmd` does not contain the saved candidate motifs. This can be computationally expensive, especially when ProbKMA-FMD is executed sequentially (`worker_number = 1`). The directory `sim300_fmd` also contains the final results of ProbKMA-FMD.

- `simulated_300_alpha05_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated in the file `simulated_length_300.R`, using partially random initialisation and alpha=0.5 (see Section "Additional details on simulation study"). The directory `sim300_alpha05_fmd` contains the candidate motifs resulting from ProbKMA-FMD. This allows to reduce the computational cost by directly uploading the candidate motifs without re-running all steps of ProbKMA-FMD. The code performs a new run of ProbKMA-FMD in the case the result directory `sim300_alpha05_fmd` does not contain the saved candidate motifs. This can be computationally expensive, especially when ProbKMA-FMD is executed sequentially (`worker_number = 1`). The directory `sim300_alpha05_fmd` also contains the final results of ProbKMA-FMD.

- `simulated_300_alpha01_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated in the file `simulated_length_300.R`, using partially random initialisation and alpha=0.1 (see Section "Additional details on simulation study"). The directory `sim300_alpha01_fmd` contains the candidate motifs resulting from ProbKMA-FMD. This allows to reduce the computational cost by directly uploading the candidate motifs without re-running all steps of ProbKMA-FMD. The code performs a new run of ProbKMA-FMD in the case the result directory `sim300_alpha01_fmd` does not contain the saved candidate motifs. This can be computationally expensive, especially when ProbKMA-FMD is executed sequentially (`worker_number = 1`). The directory `sim300_alpha01_fmd` also contains the final results of ProbKMA-FMD.

- `simulated_300_only_random_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated in the file `simulated_length_300.R`, using only random initialisation and alpha=0.9 (see Section "Additional details on simulation study"). The directory `sim300_only_random_fmd` contains the candidate motifs resulting from ProbKMA-FMD. This allows to reduce the computational cost by directly uploading the candidate motifs without re-running all steps of ProbKMA-FMD. The code performs a new run of ProbKMA-FMD in the case the result directory `sim300_only_random_fmd` does not contain the saved candidate motifs. This can be computationally expensive, especially when ProbKMA-FMD is executed sequentially (`worker_number = 1`). The directory `sim300_only_random_fmd` also contains the final results of ProbKMA-FMD.

## Empirical analysis

The folder `Empirical_analysis`  encloses the following files and folders:

- `SP500_wrds.xlsx ` - the file with daily prices of the considered stocks downloaded from Wharton Research Data Services platform.

- `pre-processing.R` - R code for detrending and implementing the hierarchical clustering using the Spearman coefficient on the same time periods in order to select the less correlated input curves (the less correlated stocks) for probabilistic K-means (probKMA). Afterwards those curves are detrended and smoothed throughout their entire length to be used for functional motif discovery with ProbKMA. We save them in `data_wrds.Rdata` file. This R code also prepares the dataset that has been used for Functional Motif Search on all the curves, after the FMD has been run on the less correlated curves. It is illustrated how we detrend and smooth the daily market prices of 51 considered stocks and save them in a `wrds_all_data.Rdata` file.

- `fmd_all.R` - R code for FMD and search using the random and the partially random initializations. The datasets `data_wrds.Rdata` and `wrds_all_data.Rdata` need to be loaded, therefore the file `pre-processing.R` should be executed beforehand in order to create them. In this case in order to reduce the computational time, the directory needs to be set in the folder `stock_prices_motifs` which contains the candidate motifs resulting from FMD using random and partially random initializations. The code will perform a new FMD in the case the working directory does not contain the saved candidate motifs. In this case the process of FMD can be particularly computationally expensive when executed sequentially.
In order to obtain the forecastings using different overlaps, the parameter `overlap` in the function `motifs_forecast` needs to be changed (the value we use in the study is of 40%).

- `stock_prices_motifs` contains the candidate motifs obtained using the random and the partially random initializations when implementing the functional motif discovery with the minimum lengths of the clusters of 40,50 and 60. 

- `fmd_all_shorter.R` - R code for FMD and search using the random and the partially random initializations with the minimum motif lengths of 20, 30, 40, 50 and 60 days. In this case in order to reduce the computational time, the directory needs to be set in the folder `stock_prices_motifs_shorter` which contains the candidate motifs resulting from FMD using random and partially random initializations. The code will perform a new FMD in the case the working directory does not contain the saved candidate motifs. In this case the process of FMD can be particularly computationally expensive when executed sequentially.

- `stock_prices_motifs_shorter` contains the candidate motifs obtained using the random and the partially random initializations when implementing the functional motif discovery with the minimum lengths of the clusters of 20,30,40,50 and 60. 





## Extended probKMA-FMD

The file `fmd_functions.R` contains all R functions of our extended probKMA-FMD for discovering motifs in stock market prices (for the original version of probKMA-FMD, see [this GitHub page](https://github.com/marziacremona/ProbKMA-FMD)). In particular:

- `probKMA`: probabilistic k-means with local alignment to find candidate motifs.
      
- `probKMA_plot`: plot the results of probKMA.
      
- `find_candidate_motifs`: run multiple times probKMA function with different K, c and initializations, with the aim to find a set of candidate motifs.
      
- `filter_candidate_motifs`: filter the candidate motifs on the basis of a threshold on the average silhouette index and a threshold on the size of the curves in the motif.
      
- `cluster_candidate_motifs`: determine a global radius, group candidate motifs based on their distance, and determine a group-specific radius.
      
- `cluster_candidate_motifs_plot`: plot the results of cluster_candidate_motifs.
      
- `motifs_search`: find occurrences of the candidate motifs in the curves and sort them according to their frequencies and radius.
      
- `motifs_search_plot`: plot the results of motifs_search.

- `motifs_search_plot_norm_time`: plot the results of normalized motifs with the actual dates on the x-axis.

- `motifs_init`: function for initializing the motifs in the form of normal triangles.

- `motifs_init_rev`: function for initializing the motifs in the form of reversed triangles.

- `motifs_line`: function for initializing the motifs in the form of increasing line.

- `motifs_line_rev`: function for initializing the motifs in the form of decreasing line.

