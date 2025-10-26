
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(testthat)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
a1 <- 10
lambda1 <- 2

a2 <- -10
lambda2 <- 1

test_that("soft_c works", {
  expect_equal(soft_c(a1, lambda1), soft(a1, lambda1))    
  expect_equal(soft_c(a2, lambda2), soft(a2, lambda2))
})


# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
test_that("lasso_c works", {
  X1 <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  Y1 <- c(1, 2, 3)
  beta1 <- c(0.5, 0.25)
  lambda1 <- 0.1
  
  X2 <- matrix(c(4, 2, 6, 9, 15, 6), ncol = 2)
  Y2 <- c(4, 5, 1)
  beta2 <- c(0.75, 0.25)
  lambda2 <- 0.2
  
  expect_equal(lasso(X1, Y1, beta1, lambda1), lasso_c(X1, Y1, beta1, lambda1))
  expect_equal(lasso(X2, Y2, beta2, lambda2), lasso_c(X2, Y2, beta2, lambda2))
})

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
