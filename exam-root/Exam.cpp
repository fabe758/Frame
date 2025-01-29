
/*
 *  Exam.cpp
 *
 *
 *  Created by F. Abe
 *
 */
// include headers
// #include "../src/Lcurve.h"
#include "../src/RT-Graph.h"

// Parameter for TCanvas
RTCanvas rtc;

// Definition of quad lenses
dMlens ml41 = {dLens(0.9, dV2d(0.0, 0.0)), dLens(0.04, dV2d(1.13, 0.11)),
               dLens(0.04, dV2d(0.98, -0.21)), dLens(0.02, dV2d(1.22, -0.22))};
dMlens ml42 = {dLens(0.993, dV2d(0.0, 0.0)), dLens(0.003, dV2d(1.13, 0.21)),
               dLens(0.003, dV2d(-0.80, -0.21)),
               dLens(0.001, dV2d(1.15, -0.32))};

// Diffinition of sources
dSource src01(dCircle(dV2d(0.1, 0.00), 0.01));
dSource src001(dCircle(dV2d(0.1, 0.00), 0.001));

// Definition of FrTri
FrTri<double> fr41(ml41, dSource(dCircle(dV2d(0.1, 0.00), 0.01)));

// Definition of t0
Time<double> t0(0.0);

// Definition of motion
dLinear ln1(-0.3, t0, M_PI / 2.0, 20.0);

// definition of source motion
dSrcMotion sm01(src01, ln1);

// definition of light curve
dLcurve lc41_01({sm01}, dMlMotion(ml41), t0.range(200, -20.0, 20.0));
