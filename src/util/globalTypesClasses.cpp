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

void ElementProperties::Initialize_ElementProperties()
{
    c = sqrt(E / rho);
    Z = c * rho;
    time_e = hE / c;
}
