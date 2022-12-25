#include "globalFunctions.h"
#include <float.h>

int SetNewtonCotes_Points_AndWeights(int numSpatialSubsegmentsPerSegment,
                                     vector<double> &spatialIntegrationWeights,
                                     vector<double> &spatialIntegrationPoints) {
  int numSpatialPointsPerSegment = numSpatialSubsegmentsPerSegment + 1;
  spatialIntegrationWeights.resize(numSpatialPointsPerSegment);

  spatialIntegrationPoints.resize(numSpatialPointsPerSegment);
  double h = 1.0 / numSpatialSubsegmentsPerSegment;
  for (int i = 0; i < numSpatialPointsPerSegment; ++i)
    spatialIntegrationPoints[i] = i * h;

  if (numSpatialPointsPerSegment == 2) {
    spatialIntegrationWeights[0] = 0.5;
    spatialIntegrationWeights[1] = 0.5;
  } else if (numSpatialPointsPerSegment == 3) {
    spatialIntegrationWeights[0] = 1.0 / 6.0;
    spatialIntegrationWeights[1] = 4.0 / 6.0;
    spatialIntegrationWeights[2] = 1.0 / 6.0;
  } else if (numSpatialPointsPerSegment == 4) {
    double tmp = 1.0 / 8.0;
    spatialIntegrationWeights[0] = tmp;
    spatialIntegrationWeights[1] = 3.0 * tmp;
    spatialIntegrationWeights[2] = 3.0 * tmp;
    spatialIntegrationWeights[3] = tmp;
  } else if (numSpatialPointsPerSegment == 5) {
    double tmp = 1.0 / 90.0;
    spatialIntegrationWeights[0] = 7.0 * tmp;
    spatialIntegrationWeights[4] = 7.0 * tmp;
    spatialIntegrationWeights[1] = 32.0 * tmp;
    spatialIntegrationWeights[3] = 32.0 * tmp;
    spatialIntegrationWeights[2] = 12.0 * tmp;
  } else if (numSpatialPointsPerSegment == 6) {
    double tmp = 1.0 / 288.0;
    spatialIntegrationWeights[0] = 19.0 * tmp;
    spatialIntegrationWeights[5] = 19.0 * tmp;
    spatialIntegrationWeights[1] = 75.0 * tmp;
    spatialIntegrationWeights[4] = 75.0 * tmp;
    spatialIntegrationWeights[2] = 50.0 * tmp;
    spatialIntegrationWeights[3] = 50.0 * tmp;
  } else if (numSpatialPointsPerSegment == 6) {
    double tmp = 1.0 / 840.0;
    spatialIntegrationWeights[0] = 41.0 * tmp;
    spatialIntegrationWeights[6] = 41.0 * tmp;
    spatialIntegrationWeights[1] = 216.0 * tmp;
    spatialIntegrationWeights[5] = 216.0 * tmp;
    spatialIntegrationWeights[2] = 27.0 * tmp;
    spatialIntegrationWeights[4] = 27.0 * tmp;
    spatialIntegrationWeights[3] = 272.0 * tmp;
  } else {
    double ovsh = 1.0 / 6.0 * h;
    spatialIntegrationWeights[0] = ovsh;
    int numPtsHalf = (numSpatialPointsPerSegment - 1) / 2;
    int cntr = 1;
    for (int i = 0; i < numPtsHalf; ++i) {
      spatialIntegrationWeights[cntr++] = 4.0 * ovsh;
      spatialIntegrationWeights[cntr++] = 2.0 * ovsh;
    }
    if (numSpatialPointsPerSegment % 2 == 0) // even segments
    {
      spatialIntegrationWeights[numSpatialSubsegmentsPerSegment] = 0.5 * h;
      spatialIntegrationWeights[numSpatialSubsegmentsPerSegment - 1] =
          ovsh + 0.5 * h;
    } else
      spatialIntegrationWeights[numSpatialSubsegmentsPerSegment] = ovsh;

    double sm = 0.0;
    for (int i = 0; i < numSpatialPointsPerSegment; ++i)
      sm += spatialIntegrationWeights[i];
    if (fabs(sm - 1.0) > 1e-5) {
      for (int i = 0; i < numSpatialPointsPerSegment; ++i)
        cout << spatialIntegrationWeights[i] << '\t';
      cout << "\nsum\t" << sm << '\n';
      THROW("Sum is not one\t");
    }
  }
  return numSpatialPointsPerSegment;
}

double computeRatio(double numerator, double denominator) {
  double max_ret = 1e40;
  if (fabs(denominator) < 2 * DBL_MIN)
    if (fabs(numerator) < 2 * DBL_MIN)
      return 0;
    else if (((numerator > 0) && (denominator > 0)) ||
             ((numerator < 0) && (denominator < 0)))
      return max_ret;
    else
      return -max_ret;

  double ratio = numerator / denominator;

  if (ratio > max_ret)
    return max_ret;
  if (ratio < -max_ret)
    return -max_ret;

  return ratio;
}
