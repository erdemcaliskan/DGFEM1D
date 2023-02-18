#ifndef INHOMOGENEOUS_FIELD__H
#define INHOMOGENEOUS_FIELD__H

#include "RandomVariable.h"

using namespace std;

// one inhomogeneous field for 1D domain
class OneIHField
{
public:
	OneIHField();
	~OneIHField();

	// inData: input for data
	// inConfig: read INPUTs below (all booleans given) plus bool num_Vals_and_x_Provided, bool containsRepeatingEndPeriodicVal
	// the last 3 values are provided as pointers to provide (isPeriodic, xM, xm) in case they are provided from outside
	void Read_Initialize_OneIHField(istream& inData, istream* inConfigPtr, bool* isPeriodicPtr = NULL, double* xMPtr = NULL, double* xmPtr = NULL, int resolutionFactor = 1, setStatOp_type sso = sso_mean_arithmetic);
	// resolutionFactor 
	//					== 0 or +/-1 nothing happens
	//					>  1  -> number of segments is DECREASED by this factor (e.g. if resolutionFactor ==  10 and numSegments = 1000 -> numSegments becoms 100)
	//					<  -1 -> number of segments is INCREASED by this factor (e.g. if resolutionFactor == -10 and numSegments = 1000 -> numSegments becoms 10000)
	void Modify_Resolution(int resolutionFactor, setStatOp_type sso);
	inline bool IsEmpty() const {	return (numSegments == 0);	}

	// this gives vals[index]
	double getVertexValueByIndex(unsigned int index) const;
	double getSegmentValueByIndex(unsigned int index) const;
	double getValueByCoord(double x) const;
	unsigned int getNumVertices() const { return numVertices; }
	unsigned int getNumSegments() const { return numSegments; }
	unsigned int getNumValues() const { return vals.size(); }
	void get_domain_range(double& xm_out, double& xM_out) const { xm_out = xm; xM_out = xM; }
	void get_xs(vector<double>& xsOut) const { xsOut = xs; }
	double get_xm() const { return xm; }
	double get_xM() const { return xM; }

	void Output(ostream& out);
	OneIHField(const OneIHField& other);
	OneIHField& operator=(const OneIHField& other);

private:
	bool valsAtVertices;
	// for uniform grid, only values are read from the file and xs are uniformly distributed from xm, xM
	// if the grid is not uniform, xs are also read from the file
	bool uniformGrid;
	bool isPeriodic;
	// domain minimum and maximum x
	double xm, xM;
	gRandVar* randVariableType;

	/////////////////////////////////////////////////////////////////////////////////
	// read from a file
	// values are given at vertices or are constant in elements connecting vertices based on valsAtVertices
	vector<double> vals;


	/////////////////////////////////////////////////////////////////////////////////
	// values finally set
	// vertex positions
	vector<double> xs;
	unsigned int numVertices;
	unsigned int numSegments;

private:
	bool num_Vals_and_x_Provided;
	bool containsRepeatingEndPeriodicVal;

	void Read_InstructionsOnly(istream* inConfigPtr, bool* isPeriodicPtr = NULL, double* xMPtr = NULL, double* xmPtr = NULL);
	void Read_DataOnly(istream& inData);

	// only reads vals and xs if uniformGrid == false
	// num_Vals_Provided: the size of number of values to be read given
	// containsRepeatingEndPeriodicVal: for periodic data if values are at vertices, the input file contains the repeating value at the end
	void Read_Vals_xs(istream& in, bool num_Vals_and_x_Provided, bool containsRepeatingEndPeriodicVal);
	void Finalize_spatialPositions();
	// is randomVariableType != NULL, the assumption is that read values are standard normal and they are mapped to target random type
	void Finalize_Values();
};

void TestInhomogeneousField(string baseNameWOExt = "TestFiles/file_v1_sz0_1");

#endif