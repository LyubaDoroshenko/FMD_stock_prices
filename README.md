# Functional Motif Discovery for Stock Market Prices
This repository contains the R code and data files used in the article [Functional motif discovery in stock market prices](https://dx.doi.org/10.2139/ssrn.4642040) by Marzia A. Cremona, Lyubov Doroshenko and  Federico Severino.

The folders `Simulation_study` and `Empirical_analysis` enclose codes of the simulation study and empirical analysis, respectively. 

## Simulation study

The folder `Simulation study` encloses the following files and folders:

-`simulated_length_300.R` - R code for simulating 20 curves of length 300 following the Mixed Autoregressive MAR(1,1) process described in the Section "Simulation study". The resulting simulated data are in the file `data_smoothed.RData`.

-`simulated_300_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated in the file `simulated_length_300.R`, using partially random initialisation and alpha=0.9 (see Section "Simulation study"). The directory `sim300_fmd` contains the candidate motifs resulting from ProbKMA-FMD. This allows to reduce the computational cost by directly uploading the candidate motifs without re-running all steps of ProbKMA-FMD. The code performs a new run of ProbKMA-FMD in the case the result directory `sim300_fmd` does not contain the saved candidate motifs. This can be computationally expensive, especially when ProbKMA-FMD is executed sequentially (`worker_number = 1`). The directory `sim300_fmd` also contains the final results of ProbKMA-FMD.

-`simulated_300_alpha05_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated in the file `simulated_length_300.R`, using partially random initialisation and alpha=0.5 (see Section "Additional details on simulation study"). The directory `sim300_alpha05_fmd` contains the candidate motifs resulting from ProbKMA-FMD. This allows to reduce the computational cost by directly uploading the candidate motifs without re-running all steps of ProbKMA-FMD. The code performs a new run of ProbKMA-FMD in the case the result directory `sim300_alpha05_fmd` does not contain the saved candidate motifs. This can be computationally expensive, especially when ProbKMA-FMD is executed sequentially (`worker_number = 1`). The directory `sim300_alpha05_fmd` also contains the final results of ProbKMA-FMD.

-`simulated_300_alpha01_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated in the file `simulated_length_300.R`, using partially random initialisation and alpha=0.1 (see Section "Additional details on simulation study"). The directory `sim300_alpha01_fmd` contains the candidate motifs resulting from ProbKMA-FMD. This allows to reduce the computational cost by directly uploading the candidate motifs without re-running all steps of ProbKMA-FMD. The code performs a new run of ProbKMA-FMD in the case the result directory `sim300_alpha01_fmd` does not contain the saved candidate motifs. This can be computationally expensive, especially when ProbKMA-FMD is executed sequentially (`worker_number = 1`). The directory `sim300_alpha01_fmd` also contains the final results of ProbKMA-FMD.

## Empirical analysis

The folder `Empirical_analysis`  encloses the following files and folders:

## Extended probKMA-FMD

The file `fmd_functions.R` contains all R functions of our extended probKMA-FMD. In particular:

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

