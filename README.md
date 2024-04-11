# methylscore 
Calculation of epigenetic scores using genome-wide methylation data

## DESCRIPTION
The process requires two input files:

### Methylation data of study sample
The input should have the following format:
  - A column "cpg" containing all CpG site names (e.g., "cg04309214").
  - Individual columns with QCed methylation beta values for each site, labeled with individual IDs.
A sample file, betamatrix_test.txt, is included in the test directory.

### EWAS summary statistics
The input should have the following format:
  - A column "cpg" containing all CpG site names (e.g., "cg04309214").
  - A column "beta" containing CpG effect sizes.
  - A column "p" containing associated p-values.
A sample file, sumstats_test.txt, is included in the test directory.

## Contact
If you have any questions, please contact agonzalezsegura@ub.edu
