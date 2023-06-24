#include "globalTypesClasses.h"
#include "commonMacros.h"

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
