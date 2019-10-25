# Lensing corrections on galaxy-lensing cross correlations and galaxy-galaxy correlations 

Scripts and notebooks to compute higher order lensing corrections following http://arxiv.org/abs/1910.06722

## Description
One usually assumes that lensed fields can be modeled in the Born approximation, i.e. that one can ignore multiple deflections and simply add all lenses along the line of sight to get an accurate estimate of the lensing effect on an observed object.

![postborn_plot](/plots/lensing_illu.png)


In this repository, we provide code to estimate the size of the lowest order correction terms to this approximation focusing on two observables:
- cross correlations between galaxies and lensing (for arbitrary galaxy redshift distributions and lensing source redshifts)
- correlations between galaxy samples (for arbitrary redshift distributions)

If you are an observer and want to check if you can safely ignore these higher order corrections for your data analysis, you can use the provided notebooks to estimate their size for your specific configuration.

## Usage
To compute corrections to cross correlations between lensing convergence (at arbitrary source redshifts) and galaxy counts (LSST-like or Gaussian redhsift distribution) run and follow instructions in
- CrossCorrelationResults.ipynb

To compute corrections for correlations between galaxy counts (same or different redshift distribution) run and follow instructions in
- AutoCorrelationResults.ipynb

An example output:  

![lensing_corrections](/plots/gaussgal_z10_sigma4_cmblens_simple_bias_zmax10886/4paper/postbornterms_gaussgal_z10_sigma4_cmblens_simple_bias_zmax10886.png)

The computation of correlations (Cls) between quantities at different redshift (without Limber approximation) is based on the FFTlog formalism ([Assassi et al. 2017](https://arxiv.org/pdf/1705.05022.pdf)). We provide the tabulated outputs of the required decomposition in the data folder (for a fiducial cosmology). The parameter values can be found in cosmo_dict.json.


## Contributors
Chirag Modi, Emanuele Castorina

## Requirements
- Python>=3
- Numpy
- Scipy
- Mpi4py  

for creating a conda environment:  
```
conda env create -f postborn_conda.yml
```

## Comments
are welcome :)

