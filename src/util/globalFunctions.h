#ifndef GLOBAL_FUNCTIONS__H
#define GLOBAL_FUNCTIONS__H

#include "globalMacros.h"

// return is numSpatialPointsPerSegment
// points are between 0 and 1

int SetNewtonCotes_Points_AndWeights(int numSpatialSubsegmentsPerSegment,
                                     vector<double> &spatialIntegrationWeights,
                                     vector<double> &spatialIntegrationPoints);

double computeRatio(double numerator, double denominator);

#endif