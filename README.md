# MSFAST Bayesian (Multivariate) FPCA

Supporting materials for MSFAST: Multivariate, Sparse Bayesian FPCA using the FAST Approach

## Description

Includes STAN code, R scripts performing comparative simulations for performance evaluation, and a simple vignette illustrating application of MSFAST to a real data set of sparsely-observed child growth measures (CONTENT). The simulation scripts are divided into univariate (P = 1) and multivariate (P = 3) scenarios, as there are different comparator implementations available for the two. The vignette includes exploratory data visualization, organization of the data into the appropriate input format, fitting the STAN model, assessing convergence and global variance explained, visualizing the FPC/eigenvalue estimates with corresponding credible intervals, and performing efficient dynamic prediction.

## Getting Started

### R Package Dependencies for Comparator Methods

* [mfaces](https://cran.r-project.org/package=mfaces)
* [MFPCA](https://cran.r-project.org/package=MFPCA)
* [bayesFPCA](https://github.com/hruffieux/bayesFPCA)
* [face](https://cran.r-project.org/package=face)
* [fdapace](https://cran.r-project.org/package=fdapace)

### R Package Dependencies for Orthogonal Spline Bases

* [Splinets](https://cran.r-project.org/package=Splinets)
* [orthogonalsplinebasis](https://cran.r-project.org/package=orthogonalsplinebasis)
* [orthopolynom](https://cran.r-project.org/package=orthopolynom)

### Installing

* Clone all files from this Repository
* Ensure all general package dependencies are satisfied (Libs.R under Gen_Funcs)

### Executing program

* All simulation scripts expect to be called from Command Line. Find example below, where 10 iterations of the multivariate simulation are run with each covariate having between 5 and 15 observations:
```
Rscript Multivar_Accuracy.R "Multi_5to15" "output" 10 5 15
```
* Visualization scripts will require that simulations have been run

## Authors

[Joseph Sartini](https://jsartini.github.io/Sartini-Stats/)

## Version History

* 0.1
    * Initial Release
