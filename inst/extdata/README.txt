yll-example.rds is a complete estimate_yll() result generated from all 5,000
synthetic participants in the packaged yll_toy data. It contains no observed
participant records.

The target population is unexposed to hypertension (reference No, exposed Yes).
Covariates are sex, education_years, bmi and smoke_binary. Starting ages are
40 through 90 years in five-year steps, with an upper age of 90. The result uses
1,000 participant bootstrap replicates, seed 20260917, normal 95% pointwise
confidence intervals, and the two-worker future multisession backend.

This file lets documentation builds draw confidence intervals without repeating
the bootstrap. Read it with:
readRDS(system.file("extdata", "yll-example.rds", package = "estimandYLL"))

Reproduction and validation scripts:
https://github.com/ShuntaroS/estimandYLL/tree/main/notes/published-example
