
# ProjectSVR

<img src="man/figures/ProjectSVR-logo.png" height="200" align="right" />

[![](https://img.shields.io/badge/devel%20version-0.1.0-green.svg)](https://github.com/JarningGau/ProjectSVR)

`ProjectSVR` is a machine learning-based algorithm for mapping the query
cells onto well-constructed reference atlas.

<img src="man/figures/ProjectSVR-workflow.png" width="600" />

## Installation

Install the development version from GitHub use:

``` r
install.packages("devtools")
devtools::install_github("JarningGau/ProjectSVR")
```

## Usage/Demos

### Tutorials

-   Mapping PBMC dataset onto pre-build PBMC reference (from DISCO
    database).

-   Mapping tumor-infiltrated T cell to pan-cancer T cell landscape.

### Downloading pre-built references:

-   You can download pre-built references from [Zenodo]().

## Installation notes

`ProjectSVR` has been successfully installed on Linux using the devtools
package to install from GitHub.

Dependencies:

-   R &gt;= 4.1
-   reticulate
-   AUCell

## Reproducing results from manuscript

Code to reproduce results from the Gao et al. manuscript is available on
github.com/jarninggau/ProjectSVR\_reproducibility.
