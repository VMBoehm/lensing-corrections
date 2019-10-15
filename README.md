# Lensing corrections on galaxy-lensing cross correlations and galaxy-galaxy auto correlations 

## Description
Scripts and notebooks to compute higher order lensing corrections following (todo: link to publication)

To compute corrections to cross correlations between lensing convergence (at arbitrary source redshifts) and galaxy counts (LSST-like or Gaussian redhsift distribution) run and follow instructions in
- CrossCorrelationResults.ipynb

To compute corrections for correlations between galaxy counts (same or different redshift distribution) run and follow instructions in
- AutoCorrelationResults.ipynb

The computation of correlations (Cls) between quantities at different redshift (without Limber approximation) is based on the FFTlog formalism ([Assassi et al. 2017](https://arxiv.org/pdf/1705.05022.pdf)). We provide the tabulated outputs of the required decomposition in the data folder (for a fiducial cosmology). The parameter values can be found in cosmo_dict.json.


## Contributors
Chirag Modi(@modichirag), Emanuele Castorina

## Requirements
- Python>=3
- Numpy
- Scipy

