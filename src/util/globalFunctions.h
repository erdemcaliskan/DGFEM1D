#ifndef GLOBAL_FUNCTIONS_NEW_H
#define GLOBAL_FUNCTIONS_NEW_H

#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

// return is numSpatialPointsPerSegment
// points are between 0 and 1

int SetNewtonCotes_Points_AndWeights(int numSpatialSubsegmentsPerSegment, vector<double> &spatialIntegrationWeights,
                                     vector<double> &spatialIntegrationPoints);

double computeRatio(double numerator, double denominator);

void ReadVectorDouble(istream& in, vector<double>& dat);

bool DoublesAreEqual(double d1, double d2, double tol);

#endif