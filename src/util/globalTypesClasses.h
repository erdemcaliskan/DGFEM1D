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

typedef enum
{
    cfem_1D,
    DG_2FUV,
    DG_1F_vStar,
    DG_1F_uStar,
    WeakFormulationType_SIZE
} WeakFormulationType;

typedef enum
{
    dg_eps_m1 = -1,
    dg_eps_0,
    dg_eps_p1,
    DG_epsilonT_SIZE
} DG_epsilonT;

#endif
