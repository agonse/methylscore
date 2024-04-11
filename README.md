# methylscore 
Calculation of epigenetic scores using genome-wide methylation data


## DESCRIPTION
### Methylation data of study sample
The input containing the methylation information of the study sample requires the following format:
  - A column "cpg": contains all cpg site names, e.g. "cg04309214"
  - A column for each individual reporting the QCed methylation beta for each site. Add ID code as column name.
A dummy betamatrix_test.txt is provided in the test directory


### EWAS summary statistics
The input containing the size effect of CpGs retrieved from EWAS requires the following format:
  - A column "cpg": contains all cpg site names, e.g. "cg04309214"
  - A column "beta": contais CpG effect sizes
  - A column "p": contains the associated p-value
A dummy sumstats_test.txt is provided in the test directory

## Contact
Should any questions arise, you can contact agonzalezsegura@ub.edu.
