#ifndef GLOBAL_TYPE_CLASSES__H
#define GLOBAL_TYPE_CLASSES__H

#include "globalMacros.h"

// characteristics is the case that incoming characteristic (zero for transmitting) from outside is given
// bct_Symmetric: for a given direction, the characteristics are equal from opposite side
// bct_AntiSymmetric: for a given direction, the characteristics are opposites from opposite side
//		For example in 2D for mode I , we need bct_Symmetric for dir 0 and bct_AntiSymmetric for dir 1 (plus dir 2 for
//3D) 						  for mode II, we need bct_AntiSymmetric for dir 0 and bct_Symmetric for dir 1
typedef enum
{
    bct_Dirichlet,
    bct_Neumann,
    bct_Characteristics,
    bct_Symmetric,
    bct_AntiSymmetric,
    bct_PeriodicOrBloch,
    bct_Unspecified,
    BoundaryConditionT_SIZE
} BoundaryConditionT;

string getName(BoundaryConditionT dat);
void name2Type(string &name, BoundaryConditionT &typeVal);
ostream &operator<<(ostream &out, BoundaryConditionT dat);
istream &operator>>(istream &in, BoundaryConditionT &dat);

typedef enum {mpt_uniform, mpt_random_field, mpt_layered, matPropT_SIZE} matPropT;

string getName(matPropT dat);
void name2Type(string &name, matPropT &typeVal);
ostream &operator<<(ostream &out, matPropT dat);
istream &operator>>(istream &in, matPropT &dat);


class ElementProperties
{
  public:
    // mechanical properties
    // E, rho, damping, he, and damping given, the other computed
    void Initialize_ElementProperties();
    // inputs
    double E, rho, damping;
    double hE;

    // computed
    // c = sqrt(E/rho)
    // Z = c * rho
    // time_e = hE / c
    double c, Z, time_e;
};

typedef enum
{
    so_Riemann,
    so_Central,
    so_Alternating_sL,
    so_Alternating_sR,
    StarOption_SIZE
} StarOption;

string getName(StarOption dat);
void name2Type(string &name, StarOption &typeVal);
ostream &operator<<(ostream &out, StarOption dat);
istream &operator>>(istream &in, StarOption &dat);

typedef enum
{
    cfem,
    DG_2FUV,
    DG_1F_vStar,
    DG_1F_uStar,
    WeakFormulationType_SIZE
} WeakFormulationType;

string getName(WeakFormulationType dat);
void name2Type(string &name, WeakFormulationType &typeVal);
ostream &operator<<(ostream &out, WeakFormulationType dat);
istream &operator>>(istream &in, WeakFormulationType &dat);

typedef enum
{
	// matrix only
	smt_Mats_Only,	 // M, C, K computed
	smt_Mats_OnlyNaturalMode, // (-M omega^2 + i omega C + K) a  = 0 is solved for omega
	smt_Mats_OnlyBloch,		 // option for Bloch mode analysis 
	// matrix and RHS F
	smt_MatsF_Static, // Ka = F
	smt_MatsF_Dynamic, // Ma'' + Ca' + Ka = F is solved dynamically - time stepping option should be provided separately
	smt_MatsF_Helmholtz, // (-M omega^2 + i omega C + K) a = F
	SolutionModeT_SIZE
} SolutionModeT;

string getName(SolutionModeT dat);
void name2Type(string &name, SolutionModeT &typeVal);
ostream &operator<<(ostream &out, SolutionModeT dat);
istream &operator>>(istream &in, SolutionModeT &dat);
void SetBooleans_SolutionMode(SolutionModeT sm, bool& needForce, bool& isDynamic);

typedef enum
{
    dg_eps_m1 = -1,
    dg_eps_0,
    dg_eps_p1,
    DG_epsilonT_SIZE
} DG_epsilonT;


// st, st  + step, ... end (end <= en) 
// en included if (enInclusive)
// either step OR num is specified
// { vl } ONE value is read
class UniformVals
{
	friend 	istream& operator>>(istream &input, UniformVals& dat);
public:
	UniformVals();
	void setVals();
	void Read(istream& in);
	double st;
	double en;
	double step;
	// num = number of SEGMENTS (NOT the points)
	int num;
	bool enInclusive;

	vector<double> vals;
};

class VecOfVals
{
	friend 	istream& operator>>(istream &input, VecOfVals& dat);
public:
	VecOfVals();
	void add_VecOfVals(double val);
	void set_VecOfVals(vector<double>& finalValsIn);
	// if delVal < 0, it's treated as the number of segments
	// return is the number of segments
	int set_UniformVecOfValues(double minVal, double maxVal, double& delVal);
	void AddVecOfVals(VecOfVals& other);
	inline bool IsEmpty() const { return (finalVals.size() == 0); }
	vector<UniformVals> uniformRanges;
	vector<double> finalVals;
	double minStep; // minimum of step between different point values
	inline double& operator[] (unsigned int i) { return finalVals[i]; }
	double operator[] (unsigned int i) const { return finalVals[i]; }
	// functions needed for finding line intersection with time slices, etc.
	// returns the first index for which the value in finalVals[index] >= val within tolerance relTol * minStep
	// assumes values in finalVals are ascending order
	int FindFirstIndex_WithEqualOrLargerVal(double val, double relTol = 1e-5);
	// indexMin is the first index for which finalVals is >= minV
	// indexMax is the last index for which finalVals is <= maxV
	// used for line intersection, ...
	// returns true if indexMin <= indexMax (i.e. there is some intersection)
	// assumes values in finalVals are ascending order
	bool FinalMinMaxIndices_Of_ValsBetweenMinMax(double minV, double maxV, int& indexMin, int& indexMax, double relTol = 1e-5);
	double GetLastValue() const;
	// aux functions
	void SetMinStep();
};

#endif
