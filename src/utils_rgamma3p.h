// This file is part of the mGHS sampler by Claudio Busatto
// (https://github.com/cbusatto/mGHS, commit 1c9fd70), licensed GPL (>= 2).
// Modified for multiGGM; every change is listed in inst/mghs/UPSTREAM.md.

/* File: utils_rgamma3p.h */

#ifndef utils_rgamma3p_H
#define utils_rgamma3p_H

// multiGGM: cap on rejection-sampler attempts. On reaching it the samplers
// return NaN, which the engine reports as a stuck chain.
#define MGHS_MAX_TRIES 1000000

// multiGGM: count of extra rejection-sampler tries, reported per fit.
extern double mghs_reject_tries;

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Parabolic Cylinder D

double rg3p_approx (const double& a, const double& b, const int& c);
double rg3p_c1 (const double& a, const double& b);

#endif

