// This file is part of the mGHS sampler by Claudio Busatto
// (https://github.com/cbusatto/mGHS, commit 1c9fd70), licensed GPL (>= 2).
// Modified for multiGGM; every change is listed in inst/mghs/UPSTREAM.md.

/* File: utils_rgamma3p.h */

#ifndef utils_rIW_H
#define utils_rIW_H

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Parabolic Cylinder D

Eigen::MatrixXd rIW (const double& df, const Eigen::MatrixXd& S, const int& p);

#endif

