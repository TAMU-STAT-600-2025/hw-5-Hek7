
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
test_that("lassostandardized_c works", {
  set.seed(43)
  
  n <- 100
  p <- 5
  X1 <- scale(matrix(rnorm(n * p), n, p))
  beta_true1 <- c(2, -1, 0, 0, 0)
  Y1 <- X1 %*% beta_true1 + rnorm(n)
  Y1 <- scale(Y1, scale = FALSE)
  
  lambda1 <- 0.1
  beta_start1 <- rep(0, p)
  
  fit_1 <- fitLASSOstandardized(X1, Y1, lambda1, beta_start1)
  fit_c1 <- fitLASSOstandardized_c(X1, Y1, lambda1, beta_start1)
  
  X2 <- scale(matrix(rnorm(n * p), n, p))
  beta_true2 <- c(1, -2, 1, 0, 0)
  Y2 <- X2 %*% beta_true2 + rnorm(n)
  Y2 <- scale(Y2, scale = FALSE)
  
  lambda2 <- 0.5
  beta_start2 <- rep(0, p)
  
  fit_2 <- fitLASSOstandardized(X2, Y2, lambda2, beta_start2)
  fit_c2 <- fitLASSOstandardized_c(X2, Y2, lambda2, beta_start2)
  
  expect_equal(as.numeric(fit_c1), as.numeric(fit_1$beta), tolerance = 1e-3)
  expect_equal(as.numeric(fit_c2), as.numeric(fit_2$beta), tolerance = 1e-3)
})


# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################
library(bench)

n <- 750
p <- 25
X <- scale(matrix(rnorm(n * p), n, p))
beta_true <- c(2, -1, 2, rep(0, p - 3))
Y <- X %*% beta_true + rnorm(n)
Y <- scale(Y, scale = FALSE)

lambda <- 0.1
beta_start <- rep(0, p)

# Run benchmark
time <- bench::mark(
  r = fitLASSOstandardized(X, Y, lambda, beta_start),
  c = fitLASSOstandardized_c(X, Y, lambda, beta_start),
  iterations = 100,
  check = FALSE
)

# Print timing in seconds
print(time[, c("expression", "median")])

# Optionally convert median times to numeric seconds
median_times_sec <- as.numeric(time$median)
names(median_times_sec) <- time$expression
print(median_times_sec)


# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################
  
test_that("fitLASSOstandardized_seq_c works", {
  set.seed(43)
  
  n = 100 
  p = 5
  X1 = scale(matrix(rnorm(n*p), n, p))
  beta_true1 = c(1.5, -1, 0, 0, 0.5) 
  Y1 = X1 %*% beta_true1 + rnorm(n)
  Y1 = scale(Y1, scale = FALSE)
  
  lambda_seq1 = seq(.5, .05, length.out = 5)
  
  fit_seq <- fitLASSOstandardized_seq(X1, Y1, lambda_seq1)
  fit_seq_c <- fitLASSOstandardized_seq_c(X1, Y1, lambda_seq1)
  
  expect_equal(dim(fit_seq_c), dim(fit_seq$beta_mat))
  expect_equal(as.numeric(fit_seq_c), as.numeric(fit_seq$beta_mat), tolerance = 1e-2)
  
  # repeat
  set.seed(43)
  
  n2 <- 200
  p2 <- 10
  X2 <- scale(matrix(rnorm(n2 * p2), n2, p2))
  beta_true2 <- c(2, -0.5, 1, 0, 0, 0, 0, 0, .2, .1)
  Y2 <- X2 %*% beta_true2 + rnorm(n2)
  Y2 <- scale(Y2, scale = FALSE)
  
  lambda_seq2 = seq(.5, .05, length.out = 5)
  
  fit_seq_2 <- fitLASSOstandardized_seq(X2, Y2, lambda_seq2)
  fit_seq_c2 <- fitLASSOstandardized_seq_c(X2, Y2, lambda_seq2)
  
  expect_equal(dim(fit_seq_c2), dim(fit_seq_2$beta_mat))
  expect_equal(as.numeric(fit_seq_c2), as.numeric(fit_seq_2$beta_mat), tolerance = 1e-2)
})

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################
set.seed(43)
n <- 500
p <- 20
X <- scale(matrix(rnorm(n * p), n, p))
beta_true <- c(2, -1, 1, rep(0, p - 3))
Y <- X %*% beta_true + rnorm(n)
Y <- scale(Y, scale = FALSE)

lambda_seq <- seq(0.5, 0.05, length.out = 8)

bench_seq <- bench::mark(
  r = fitLASSOstandardized_seq(X, Y, lambda_seq),
  c = fitLASSOstandardized_seq_c(X, Y, lambda_seq),
  iterations = 100,
  check = FALSE
)

print(bench_seq[, c("expression", "median")])

median_times_sec_seq <- as.numeric(bench_seq$median)
names(median_times_sec_seq) <- bench_seq$expression
print(median_times_sec_seq)

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression
library(microbenchmark)

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

