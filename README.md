# FMD
This repository contains the R code and data files used in the article Functional Motif Discovery in Stock Market Prices by Marzia A. Cremona, Lyubov Doroshenko and  Federico Severino.

The folder Empirical analysis and Simulation study enclose codes of the empirical analysis and the simulation study, respectively.

`Simulation study` encloses the following files and folders:

-`simulated_length_300.R` - R code for simulating 20 curves of length 300 following the Mixed Autoregressive MAR(1,1) process described in the Simulation Study.

-`simulated_300_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated with the help of the Mixed Autoregressive MAR(1,1) process using the random and the partially random initialisation and alpha=0.9. The directory `sim300_fmd` contains the candidate motifs resulting from FMD using random and partially random initializations and it is possible to reduce the computational time by setting the directory in this folder and uploading the candidate motifs. The code will perform a new FMD in the case the working directory does not contain the saved candidate motifs. In this case the process of FMD can be computationally expensive when executed sequentially.

-`simulated_300_alpha05_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated with the help of the Mixed Autoregressive MAR(1,1) process using the random and the partially random initialisation and alpha=0.5. The directory `sim300_alpha05_fmd` contains the candidate motifs resulting from FMD using random and partially random initializations and it is possible to reduce the computational time by setting the directory in this folder and uploading the candidate motifs. The code will perform a new FMD in the case the working directory does not contain the saved candidate motifs. In this case the process of FMD can be computationally expensive when executed sequentially.

-`simulated_300_alpha01_fmd.R` - R code for running functional motif discovery (FMD) and search on the 20 curves of length 300 simulated with the help of the Mixed Autoregressive MAR(1,1) process using the random and the partially random initialisation and alpha=0.1. The directory `sim300_alpha01_fmd` contains the candidate motifs resulting from FMD using random and partially random initializations and it is possible to reduce the computational time by setting the directory in this folder and uploading the candidate motifs. The code will perform a new FMD in the case the working directory does not contain the saved candidate motifs. In this case the process of FMD can be computationally expensive when executed sequentially.
