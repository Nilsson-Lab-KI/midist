# R package for MID distance calculations

This repository provides the `midist` R package for computing distances between mass isotopomer distributions (MIDs) from mass spectrometry metabolomics data. The concept of using MID distances to discover related molecules from mass spectrometry data is described in the paper "Isotope tracing-based metabolite identification for mass spectrometry metabolomics", https://doi.org/10.1101/2025.04.07.647691 

Following R package standards, all dependencies required for this software can be found in the `NAMESPACE` file. The package has been tested with R v.4.2.3 on Windows and Mac OS.

## Installation

Using `devtools`, the `midist `package can be installed with
```
install_github('Nilsson-Lab-KI/midist')
```
Installation usually takes only a minute or two.

## Demonstration data and run times

For a demonstration of package usage on real-world and simulated data, please see https://github.com/Nilsson-Lab-KI/isotope-met-id
Run time is typically dominated by MID distance matrix computation (see function `conv_reduce_all`), which can be considerable for large data sets. For our example data set with 721 metabolites, distance matrix computation takes a few minutes on a typical laptop. However, computation time grows rapidly with increasing number of metabolites.
