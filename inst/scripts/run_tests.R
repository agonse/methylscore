## Test methylscore package
library(methylscore)

# Set paths to the test data
beta_file <- system.file("extdata", "betamatrix_test.txt", package = "methylscore")
sumstats_full  <- system.file("extdata", "sumstats", "sumstats_test.txt",
                              package = "methylscore")
ewas_path  <- dirname(sumstats_full)
ewas_file  <- basename(sumstats_full)

# Example 1: calculating a single methylation profile score
episcore <- calculate_episcore(
    thrs.criteria = 0.05,
    beta.file = beta_file,
    ewas.path = ewas_path,
    ewas.file = ewas_file,
    missingness = 0.2
)

# Print the results
print(episcore)

# Example 2: calculating a multiple methylation profile score
episcore_multiple <- calculate_multiple_episcores(
    thrs.criteria = 0.05,
    beta.file = beta_file,
    ewas.list.path = ewas_path,
    missingness = 0.2
)

# Print the results
print(episcore_multiple)
