## Test methylscore package
library(methylscore)

# Set paths to the test data
beta_path <- system.file("extdata", "betamatrix_test.txt", package = "methylscore")
sumstats_path <- system.file("extdata", "sumstats", "sumstats_test.txt", package = "methylscore")

# Example 1: calculating a single episcore
episcore <- calculate_episcore(
    thrs.criteria = 0.05,
    beta = beta_path,
    ewas.file = sumstats_path,
    missingness = 0.2
)

# Print the results
print(episcore)
