#include "globalTypesClasses.h"
#include "commonMacros.h"
#include "globalFunctions.h"
#include "Dims.h"

string getName(BoundaryConditionT dat)
{
    if (dat == bct_Dirichlet)
        return "Dirichlet";
    if (dat == bct_Neumann)
        return "Neumann";
    if (dat == bct_Characteristics)
        return "Characteristics";
    if (dat == bct_Unspecified)
        return "Unspecified";
    if (dat == bct_Symmetric)
        return "Symmetric";
    if (dat == bct_AntiSymmetric)
        return "AntiSymmetric";
    if (dat == bct_PeriodicOrBloch)
        return "PeriodicOrBloch";
    cout << (int)dat << '\n';
    THROW("invalid dat");
}

void name2Type(string &name, BoundaryConditionT &typeVal)
{
    int num;
    bool success = fromString(name, num);
    if (success)
    {
        if (num >= BoundaryConditionT_SIZE)
            THROW("too large of a number\n");
        typeVal = (BoundaryConditionT)num;
        return;
    }
    // at this point we know that name is not an integer ...
    for (int i = 0; i < BoundaryConditionT_SIZE; ++i)
    {
        typeVal = (BoundaryConditionT)
            i; // casting integer to BoundaryConditionT, if we don't cast it C++ gives a compile error
        string nameI = getName(typeVal);
        if (name == nameI)
            return;
    }
    cout << "name\t" << name << '\n';
    THROW("wrong name reading BoundaryConditionT\n");
}

// operator for output
ostream &operator<<(ostream &out, BoundaryConditionT dat)
{
    string name = getName(dat);
    out << name;
    return out;
}

// operator for input
istream &operator>>(istream &in, BoundaryConditionT &dat)
{
    string name;
    string buf;
    READ_NSTRING(in, buf, name);
    name2Type(name, dat);
    return in;
}

///////////////
string getName(matPropT dat)
{
	if (dat == mpt_uniform)
		return "uniform";
	if (dat == mpt_random_field)
		return "random_field";
	if (dat == mpt_layered)
		return "layered";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string &name, matPropT &typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= matPropT_SIZE)
			THROW("too large of a number\n");
		typeVal = (matPropT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < matPropT_SIZE; ++i)
	{
		typeVal = (matPropT)
			i; // casting integer to matPropT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading matPropT\n");
}

// operator for output
ostream &operator<<(ostream &out, matPropT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

// operator for input
istream &operator>>(istream &in, matPropT &dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}


///////////////
string getName(StarOption dat)
{
	if (dat == so_Riemann)
		return "Riemann";
	if (dat == so_Central)
		return "Central";
	if (dat == so_Alternating_sL)
		return "Alternating_sL";
	if (dat == so_Alternating_sR)
		return "Alternating_sR";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string &name, StarOption &typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= StarOption_SIZE)
			THROW("too large of a number\n");
		typeVal = (StarOption)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < StarOption_SIZE; ++i)
	{
		typeVal = (StarOption)
			i; // casting integer to StarOption, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading StarOption\n");
}

// operator for output
ostream &operator<<(ostream &out, StarOption dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

// operator for input
istream &operator>>(istream &in, StarOption &dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

///////////////
string getName(WeakFormulationType dat)
{
	if (dat == cfem)
		return "cfem";
	if (dat == DG_2FUV)
		return "DG_2FUV";
	if (dat == DG_1F_vStar)
		return "DG_1F_vStar";
	if (dat == DG_1F_uStar)
		return "DG_1F_uStar";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string &name, WeakFormulationType &typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= WeakFormulationType_SIZE)
			THROW("too large of a number\n");
		typeVal = (WeakFormulationType)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < WeakFormulationType_SIZE; ++i)
	{
		typeVal = (WeakFormulationType)
			i; // casting integer to WeakFormulationType, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading WeakFormulationType\n");
}

// operator for output
ostream &operator<<(ostream &out, WeakFormulationType dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

// operator for input
istream &operator>>(istream &in, WeakFormulationType &dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}


///////////////
string getName(SolutionModeT dat)
{
	if (dat == smt_Mats_Only)
		return "Mats_Only";
	if (dat == smt_Mats_OnlyNaturalMode)
		return "Mats_OnlyNaturalMode";
	if (dat == smt_Mats_OnlyBloch)
		return "Mats_OnlyBloch";
	if (dat == smt_MatsF_Static)
		return "MatsF_Static";
	if (dat == smt_MatsF_Dynamic)
		return "MatsF_Dynamic";
	if (dat == smt_MatsF_Helmholtz)
		return "MatsF_Helmholtz";
	cout << (int)dat << '\n';
	THROW("invalid dat");
}

void name2Type(string &name, SolutionModeT &typeVal)
{
	int num;
	bool success = fromString(name, num);
	if (success)
	{
		if (num >= SolutionModeT_SIZE)
			THROW("too large of a number\n");
		typeVal = (SolutionModeT)num;
		return;
	}
	// at this point we know that name is not an integer ...
	for (int i = 0; i < SolutionModeT_SIZE; ++i)
	{
		typeVal = (SolutionModeT)i; // casting integer to SolutionModeT, if we don't cast it C++ gives a compile error
		string nameI = getName(typeVal);
		if (name == nameI)
			return;
	}
	cout << "name\t" << name << '\n';
	THROW("wrong name reading SolutionModeT\n");
}

// operator for output
ostream &operator<<(ostream &out, SolutionModeT dat)
{
	string name = getName(dat);
	out << name;
	return out;
}

// operator for input
istream &operator>>(istream &in, SolutionModeT &dat)
{
	string name;
	string buf;
	READ_NSTRING(in, buf, name);
	name2Type(name, dat);
	return in;
}

void SetBooleans_SolutionMode(SolutionModeT sm, bool & needForce, bool & isDynamic)
{
	isDynamic = false;
	needForce = false;
	if ((sm == smt_MatsF_Static) || (sm == smt_MatsF_Helmholtz))
		needForce = true;
	else if (sm == smt_MatsF_Dynamic)
	{
		needForce = true;
		isDynamic = true;
	}
}

void ElementProperties::Initialize_ElementProperties()
{
    c = sqrt(E / rho);
    Z = c * rho;
    time_e = hE / c;
}


UniformVals::UniformVals()
{
	st = 0.0;
	en = 1.0;
	step = -1.0;
	num = -1;
	enInclusive = true;
}

void UniformVals::Read(istream& in)
{
	string buf;
	READ_NSTRING(in, buf, buf);
	if (buf != "{")
		THROW("read should start with {");
	READ_NSTRING(in, buf, buf);
	while (buf != "}")
	{
		double vl;
		if (fromString(buf, vl) == true)
		{
			st = vl;
			READ_NSTRING(in, buf, buf);
			if (buf != "}")
				THROW("Number is read\n");
			return;
		}
		else if ((buf == "st") || (buf == "start") || (buf == "val"))
		{
			READ_NDOUBLE(in, buf, st);
		}
		else if ((buf == "en") || (buf == "end"))
		{
			READ_NDOUBLE(in, buf, en);
		}
		else if (buf == "step")
		{
			READ_NDOUBLE(in, buf, step);
		}
		else if ((buf == "nu") || (buf == "num"))
		{
			READ_NINTEGER(in, buf, num);
		}
		else if (buf.compare(0, 1, "in") == 0)
		{
			READ_NBOOL(in, buf, enInclusive);
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(in, buf, buf);
	}
}

void UniformVals::setVals()
{
	if (step > 0.0)
	{
		//		num = (int)ceil((en - st) / step - 1e-7);
		num = (int)floor((en - st) / step + 1e-7);
		if (num < 1)
			THROW("num must be greater than 1");
	}
	else
	{
		if (num < 0)
		{
			vals.push_back(st);
			return;
		}
		step = (en - st) / (double)num;
	}
	int numPts = num + (int)(enInclusive);
	vals.resize(numPts);
	for (int i = 0; i < numPts; ++i)
		vals[i] = st + i * step;
}

istream& operator>>(istream &input, UniformVals& dat)
{
	dat.Read(input);
	dat.setVals();
	return input;
}

istream& operator>>(istream &input, VecOfVals& dat)
{
	int num = 1;
	int tolSz = 0;
	bool bPIflag = false;
	string buf;
	READ_NSTRING(input, buf, buf);
	if (buf != "{")
		THROW("read should start with {");
	READ_NSTRING(input, buf, buf);
	while (buf != "}")
	{
		int  nm;
		if (fromString(buf, nm) == true)
		{
			num = nm;
		}
		else if ((buf == "nm") || (buf == "num"))
		{
			READ_NINTEGER(input, buf, num);
		}
		else if ((buf == "val") || (buf == "vals"))
		{
			dat.uniformRanges.resize(num);
			for (int i = 0; i < num; ++i)
			{
				input >> dat.uniformRanges[i];
				tolSz += dat.uniformRanges[i].vals.size();
			}
		}
		else if ((buf == "pi") || (buf == "PI"))
		{
			bPIflag = true;
		}
		else
		{
			cout << "buf:\t" << buf << '\n';
			THROW("invalid option\n");
		}
		READ_NSTRING(input, buf, buf);
	}
	dat.finalVals.resize(tolSz);
	int pos = 0;
	for (int i = 0; i < num; ++i)
	{
		for (int j = 0; j < dat.uniformRanges[i].vals.size(); ++j)
		{
			dat.finalVals[pos] = dat.uniformRanges[i].vals[j];
			if (bPIflag)
				dat.finalVals[pos] *= PI;
			pos++;
		}
	}
	dat.SetMinStep();
	return input;
}

//2019
VecOfVals::VecOfVals()
{
	minStep = 0.0;
}

void VecOfVals::add_VecOfVals(double val)
{
	bool found = false;
	for (int i = 0; i < finalVals.size(); ++i)
	{
		if (DoublesAreEqual(finalVals[i], val, 1e-7) == true)
		{
			found = true;
			break;
		}
	}
	if (!found)
	{
		double relTol = 1e-6;
		int sz = finalVals.size();
		if ((sz == 0) || (val > finalVals[sz - 1]))
			finalVals.push_back(val);
		else
		{
			int ind = FindFirstIndex_WithEqualOrLargerVal(val, relTol);
			finalVals.resize(sz + 1);
			for (int i = sz; i > ind; --i)
				finalVals[i] = finalVals[i - 1];
			finalVals[ind] = val;
		}
		SetMinStep();
	}
}

void VecOfVals::set_VecOfVals(vector<double>& finalValsIn)
{
	finalVals = finalValsIn;
	SetMinStep();
}

int VecOfVals::set_UniformVecOfValues(double minVal, double maxVal, double & delVal)
{
	int sz;
	if (delVal > 0)
	{
		minStep = delVal;
		sz = floor((maxVal - minVal) / delVal + 1e-6) + 1;
	}
	else
	{
		sz = (int)floor(-delVal + 1e-6);
		minStep = (maxVal - minVal) / sz;
		++sz;
		delVal = minStep;
	}
	finalVals.resize(sz);
	for (int i = 0; i < sz; ++i)
		finalVals[i] = minVal + i * minStep;
	delVal = minStep;
	return sz;
}

void VecOfVals::AddVecOfVals(VecOfVals& other)
{
	for (int i = 0; i < other.finalVals.size(); ++i)
		add_VecOfVals(other.finalVals[i]);
}

void VecOfVals::SetMinStep()
{
	int sz = finalVals.size();
	if (sz == 0)
		return;
	if (sz == 1)
	{
		minStep = abs(finalVals[0]);
		return;
	}
	minStep = abs(finalVals[1] - finalVals[0]);
	double tmp;
	for (int i = 1; i < (sz - 1); ++i)
	{
		tmp = abs(finalVals[i + 1] - finalVals[i]);
		minStep = MIN(minStep, tmp);
	}
}

int VecOfVals::FindFirstIndex_WithEqualOrLargerVal(double val, double relTol)
{
	if (val < (finalVals[0] + relTol * minStep))
		return 0;
	for (int index = 0; index < finalVals.size(); ++index)
	{
		if (finalVals[index] >= val)
			return index;
	}
	return -1; // val is larger than all content in finalVals
}

bool VecOfVals::FinalMinMaxIndices_Of_ValsBetweenMinMax(double minV, double maxV, int& indexMin, int& indexMax, double relTol)
{
	indexMin = FindFirstIndex_WithEqualOrLargerVal(minV, relTol);
	if (indexMin < 0)
		return false; // all values are greater than minV
	maxV += relTol * minStep;
	for (indexMax = indexMin; indexMax < finalVals.size(); ++indexMax)
	{
		if (finalVals[indexMax] > maxV)
		{
			--indexMax;
			return (indexMax >= indexMin);
		}
	}
	// indexMax = finalVals.size() - 1, maxVa >= maximum in the list
	--indexMax;
	return (indexMax >= indexMin);
}

double VecOfVals::GetLastValue() const
{
	int sz = finalVals.size();
	if (sz > 0)
		return finalVals[sz - 1];
	THROW("Cannot compute final value\n");
}

