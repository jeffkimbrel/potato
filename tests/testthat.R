# This file is part of the standard R package testing structure
# It runs all tests in tests/testthat/ when R CMD check is run

library(testthat)
library(potato)

test_check("potato")
