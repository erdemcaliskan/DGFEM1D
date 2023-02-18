#ifndef LAYERED_PROPERTIES__H
#define LAYERED_PROPERTIES__H

#include "globalMacros.h"
#include <map>
using namespace std;

class Bulk_Elastic_Modifier;

class oneBulk_Elastic_Prop
{
public:
	friend istream &operator>>(istream &in, oneBulk_Elastic_Prop& dat);
	double getLength() const { return length; }
	void setLength(double lengthIn) { length = lengthIn; }
	void setFinalValue(const oneBulk_Elastic_Prop& baseMat, const Bulk_Elastic_Modifier* bemPtr);
	oneBulk_Elastic_Prop();
	double E, rho, damping;
private:
	double length;
};

class Bulk_Elastic_Modifier
{
public:
	Bulk_Elastic_Modifier();
	// returns true if it reads one instance, false if not (i.e. reading })
	bool Read_Bulk_Elastic_Modifier(istream& in);
	void Print(ostream& out) const;

	// read data
	double length; // length of the segment
	GID	bulk_flag;
	// property factors
	double CFactor, rhoFactor, dampingFactor;

	// computed data
	// if any of the factors != 1, this modifier modifies reference elastic property
	bool b_modifies;

	//auxiliary functions
	void Initialize_Bulk_Elastic_Modifier();
};

class Bulk_Elastic_Prop
{
public:
	Bulk_Elastic_Prop();
	void Read_Bulk_Elastic_Prop(istream& in, int serialNumber = -1);

	// read data
	map<GID, oneBulk_Elastic_Prop> baseBulkProperties;
	vector<Bulk_Elastic_Modifier> bulkModifiers;
	// directSpaceSizeModifier == 0 do nothing
	// > 0	-> bulk sizes should be smaller than this
	// < 0  -> refinement ratio (e.g. -2 -> each segment is divided into two)
	double directSpaceSizeModifier;
	unsigned int numRepeatSequence;

	/// computed
	vector<oneBulk_Elastic_Prop> finalBulkSegments;
	unsigned int sz_oneSequence, sz_allsequences;
private:
	bool b_directSpaceSizeModifier;
	void InitializeAfterRead();
};

#endif
