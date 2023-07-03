#ifndef RANDOM_VARIABLES__H
#define RANDOM_VARIABLES__H

#include "commonMacros.h"

#include <iostream> /* for input and output using >> and << operator */
#include <fstream>  /* file streams for input and output */
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>
#include <assert.h>     /* assert */
#include <math.h>     /* math functions */
#include <array>
#include <vector>

#ifndef NUMERIC_PINF
#define NUMERIC_PINF 1E6
#endif

#ifndef NUMERIC_NINF
#define NUMERIC_NINF -1E6
#endif

#ifndef THROW
#if VCPP
#define THROW(msg){{char tmpstr[255];sprintf_s(tmpstr, "In %s, line %d : %s \n", __FILE__, __LINE__, msg);	cerr<<tmpstr;getchar(); getchar(); throw(tmpstr); }}
#else
#define THROW(msg){{char tmpstr[255];sprintf(tmpstr, "In %s, line %d : %s \n", __FILE__, __LINE__, msg);	cerr<<tmpstr;getchar(); getchar(); throw(tmpstr); }}
#endif
#endif

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

typedef enum {pUniform, pTriangle, pNormal, pLogNormal, pWeibull,	pEmpirical, prob_distrib_type_SIZE} prob_distrib_type;

string getName(prob_distrib_type dat);
void name2Type(string& name, prob_distrib_type& typeVal);
ostream& operator<<(ostream& out, prob_distrib_type dat);
istream& operator>>(istream& in, prob_distrib_type& dat);

// computes min, max, mean (arithmetic, harmonic, geometric), sdiv, ... of a set of numbers
// sso_all_field_val returns the whole input vector
// sso_all_field_time_val is similar to above and is added as a kludge to help in DomainPostProcessingS2 and S3 extracting time history vector fiels.  sso_all_field_val  -> extracts spatial fields
// sso_number -> length of the vector
typedef enum { sso_none = -1, sso_valStart, sso_valEnd, sso_index_min, sso_min, sso_index_max, sso_max, sso_mean_arithmetic, sso_sdiv, sso_cov, sso_var, sso_mean_harmonic, sso_mean_geometric, sso_number, sso_all_field_val, sso_all_field_time_val, setStatOp_type_SIZE} setStatOp_type;
#define setStatOp_type_SIZE_SHORT	sso_mean_harmonic

string getName(setStatOp_type dat);
void name2Type(string& name, setStatOp_type& typeVal);
ostream& operator<<(ostream& out, setStatOp_type dat);
istream& operator>>(istream& in, setStatOp_type& dat);
double getStatvalue(const vector<double>& vals, setStatOp_type sso);
double getWeightedStatvalue(const vector<double>& vals, const vector<double>& weights, setStatOp_type sso);
double getWoWOWeightStatvalue(const vector<double>& vals, bool hasWeight, const vector<double>& weights, setStatOp_type sso);

typedef enum { vpo_none = -1, vpo_log10, vpo_loge, vpo_log2, valPOperT_SIZE } valPOperT;
string getName(valPOperT dat);
void name2Type(string& name, valPOperT& typeVal);
ostream& operator<<(ostream& out, valPOperT dat);
istream& operator>>(istream& in, valPOperT& dat);
double get_valPOperT(double valIn, valPOperT oper);

class statParamHolder {
public:
	statParamHolder();
	void Read(istream& in); // , double dd2 = -1.0);
	void WriteData(ostream& out);

	double minV, maxV;
	double mode;	// used fo triangular distribution, position of max likelihood in general

	double mean;
	double sdiv;

	double mu, sigma; // for normal and log normal distributions

	double scale; // lambda for weibull
	double shape; // m (k) for weibull

	// these booleans are used to store what values are read
	bool minMax_read, mode_read, mean_sdiv_read, mu_sigma_read;
};

class gRandVar {
public:
	/*Accessors*/
	virtual ~gRandVar() {};
	void ReadParameters(istream& in);
	// return value correspond to the given standard normal value for this random variable PDF
	double TurnStandardNormalValue2ThisRandom(double standardNormalValue) const;
	virtual double getRandomValue() const;
	virtual double getCDF(double x)  const = 0;
	virtual double getInverseCDF(double p) const = 0; //Quantile
	virtual double getPDF(double x)  const = 0;

	virtual void Finalize_Read();
	gRandVar* CreateCopy();
	prob_distrib_type	dist_type;
	statParamHolder paras;
};

class Uniform_RandVar : public gRandVar {
public:
	virtual double getCDF(double x) const;
	virtual double getInverseCDF(double p) const;
	virtual double getPDF(double x) const;
	virtual double getRandomValue() const;

	virtual void Finalize_Read();
	double span, inv_span;
};

class Triangle_RandVar : public gRandVar {
public:
	virtual double getCDF(double x) const;
	virtual double getInverseCDF(double p) const;
	virtual double getPDF(double x) const;

	virtual void Finalize_Read();
	// m: min, M: Max, md: mode
	// span = M - m
	// md_m = mode - m, M_md = Max - md, p_md = probablity of mode point = 2 span^-1, CCD_md = (mode - m) / (M - m) 
	double span, inv_span, md_m, inv_md_m, M_md, inv_M_md, p_md, cdf_md;

};

class Normal_RandVar : public gRandVar {
public:
	virtual double getCDF(double x) const;
	virtual double getInverseCDF(double p) const;
	virtual double getPDF(double x) const;
	virtual double getRandomValue() const;
};

class LogNormal_RandVar : public gRandVar {
public:
	virtual double getCDF(double x) const;
	virtual double getInverseCDF(double p) const;
	virtual double getPDF(double x) const;

	virtual void Finalize_Read();
};

class Weibull_RandVar : public gRandVar {
public:
	virtual double getCDF(double x) const;
	virtual double getInverseCDF(double p) const;
	virtual double getPDF(double x) const;
	virtual double getRandomValue() const;
};

gRandVar* gRandVar_Factory(prob_distrib_type dist_type);

gRandVar* ReadRandomVariableFromFile(istream& in);

double NormalCDFInverse(double p);

//ref: http://www.johndcook.com/blog/cpp_phi_inverse/
double RationalApproximation(double t);

int slrand();

#endif