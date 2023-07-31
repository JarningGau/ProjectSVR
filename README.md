
# ProjectSVR <img src="man/figures/ProjectSVR-logo.png" align="right" width=80px/>

[![](https://img.shields.io/badge/devel%20version-0.2.0-green.svg)](https://github.com/JarningGau/ProjectSVR)

`ProjectSVR` is a machine learning-based algorithm for mapping the query
cells onto well-constructed reference atlas.

<img src="man/figures/ProjectSVR-workflow.png" width="600" />

## Installation

Install the development version from GitHub use:

``` r
install.packages("devtools")
devtools::install_github("JarningGau/ProjectSVR")
```

## Related resources

### Reference atlas

The reference cell atlases involved in ProjectSVR paper are available at
<https://zenodo.org/record/8191576>.

### Query dataset

The query datasets involved in ProjectSVR paper are available at
<https://zenodo.org/record/8191576>.

### Pre-built reference model

You can download pre-built references from
[Zenodo](https://zenodo.org/record/8191559).

| Name                                                            | Source                                         | Version | Download                                                                     |
|-----------------------------------------------------------------|------------------------------------------------|---------|------------------------------------------------------------------------------|
| PBMC (DISCO)                                                    | <https://www.immunesinglecell.org/atlas/blood> | 0.1     | [download](https://zenodo.org/record/8191559/files/model.disco_pbmc.rds)     |
| Mouse testicular cell atlas (mTCA)                              | This paper                                     | 0.1     | [download](https://zenodo.org/record/8191559/files/model.mTCA.rds)           |
| Maternal-fetal interface atlas (Vento 2018)                     | <https://doi.org/10.1038/s41586-018-0698-6>    | 0.1     | [download](https://zenodo.org/record/8191559/files/model.Vento2018.MFI.rds)  |
| Pan cancer tumor infiltrated CD4+ T cell landscape (Zheng 2021) | <https://doi.org/10.1126/science.abe6474>      | 0.1     | [download](https://zenodo.org/record/8191559/files/model.Zheng.CD4Tcell.rds) |
| Pan cancer tumor infiltrated CD8+ T cell landscape (Zheng 2021) | <https://doi.org/10.1126/science.abe6474>      | 0.1     | [download](https://zenodo.org/record/8191559/files/model.Zheng.CD8Tcell.rds) |

## Usage/Demos

### Tutorials

-   [Quick start: mapping PBMC dataset onto pre-build PBMC
    reference.](https://jarninggau.github.io/ProjectSVR/articles/quick_start.html)

-   [Mapping tumor-infiltrated T cell to pan-cancer T cell
    landscape.](https://jarninggau.github.io/ProjectSVR/articles/ICB_breast_cancer.html)

-   [Mapping the genetic perturbed germ cells onto mouse testicular cell
    atlas.](https://jarninggau.github.io/ProjectSVR/articles/mTCA_perturbed_germ_cells.html)

## Installation notes

`ProjectSVR` has been successfully installed on Linux using the devtools
package to install from GitHub.

### Dependencies

-   R &gt;= 4.1

### External packages

Install `AUCell` or `UCell` for signature score calculation.

``` r
## install UCell
# R = 4.3
BiocManager::install("UCell") # or
# R < 4.3
remotes::install_github("carmonalab/UCell", ref="v1.3")
## install AUCell
BiocManager::install("AUCell")
```

We provided a wrapper
[`RunCNMF`](https://jarninggau.github.io/ProjectSVR/reference/RunCNMF.html)
of python pacakge [`cnmf`](https://github.com/dylkot/cNMF) for feature
selection. If you want to use it, you should install `cnmf` through
`reticulate`.

``` r
install.packages("reticulate")
reticulate::install_miniconda()
## install sceasy for single cell data format transformation.
devtools::install_github("cellgeni/sceasy")
reticulate::py_install("anndata")
## install cnmf package via reticulate
reticulate::py_install("cnmf")
```

## Reproducing results from manuscript

Code to reproduce results from the Gao et al. manuscript is available at
github.com/jarninggau/ProjectSVR\_reproducibility.

## Code of Conduct

Please note that the ProjectSVR project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
