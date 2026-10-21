// This file is part of the mGHS sampler by Claudio Busatto
// (https://github.com/cbusatto/mGHS, commit 1c9fd70), licensed GPL (>= 2).
// Modified for multiGGM; every change is listed in inst/mghs/UPSTREAM.md.

#include <RcppEigen.h>


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]

long long factorial (const int& n)
{
  long long f = 1;
  for (int i = 1; i <= n; ++i)
    f *= i;
  return f;
}