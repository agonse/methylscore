# methylscore 
This repositry contains the code to calculate specific signatures in DNA methylation. An epigenetic score can be created to estimate the effect of a particular event, such as an environmental factors, disorders, etc.

## DESCRIPTION
The process requires two input files:

### Methylation data of study sample (beta matrix)
The input should have the following format:
  - A column "cpg" containing all CpG site names (e.g., "cg04309214").
  - Individual columns with QCed methylation beta values for each site, labeled with individual IDs.
A sample file, betamatrix_test.txt, is included in the test directory.

### EWAS summary statistics
The input should have the following format:
  - A column "cpg" containing all CpG site names (e.g., "cg04309214").
  - A column "beta" containing CpG effect sizes.
  - A column "p" containing associated p-values.
A sample file, sumstats_test.txt, is included in the test directory. For multiple episcore construction,  sumstats_test1.txt is provided

## Contact
If you have any questions, please contact agonzalezsegura@ub.edu
