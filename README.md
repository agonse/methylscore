# methylscore 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14652220.svg)](https://doi.org/10.5281/zenodo.14652220)

The **Methylscore** repository contains R functions designed to calculate specific signatures in DNA methylation. These functions allow researchers to create epigenetic profile scores (MPS) that estimate the aassociation of environmental factors and health conditions on DNA methylation patterns.

## Description

The `methylscore` package processes two main types of input files:

### Methylation Data (Beta Matrix)

The methylation data should be formatted as follows:

- A column labeled `cpg` containing all CpG site names (e.g., `cg04309214`).
- Individual columns for each sample with QCed methylation beta values for each site, labeled with individual IDs.

A sample file, `betamatrix_test.txt`, is included in the `test` directory.

### EWAS Summary Statistics

The summary statistics file should be formatted as follows:

- A column labeled `cpg` containing all CpG site names (e.g., `cg04309214`).
- A column labeled `beta` containing CpG effect sizes.
- A column labeled `p` containing associated p-values.

For constructing multiple episcores, sample files such as `sumstats_test1.txt` are provided in the `test` directory.

## Installation

To install the `methylscore` package, use the following R command:

```r
# Install directly from GitHub
devtools::install_github("agonse/methylscore")
```

## Testing

After installing the `methylscore` package, you can test its functionality by running the included test script. This script uses example data files to perform basic operations provided by the package.

To run the test script, use the following commands in R:

```r
source(system.file("scripts", "run_tests.R", package = "methylscore"))
```

## Citation
If you use **methylscore** in your research, please cite it as follows:

> Segura, A. G. (2025). methylscore: epigenetic score calculation tool (v1.0). Zenodo. https://doi.org/10.5281/zenodo.14652220

By citing methylscore, you acknowledge the work and effort that went into developing this software and help others discover it.

## Contact
If you have any questions, please contact agonzalezsegura@ub.edu
