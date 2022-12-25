#ifndef GLOBAL_TYPE_CLASSES__H
#define GLOBAL_TYPE_CLASSES__H

#include "globalMacros.h"

class ElementProperties {
public:
  // mechanical properties
  // inputs
  double E, rho, damping;
  double hE;

  // computed
  // c = sqrt(E/rho)
  // Z = c * rho
  // time_e = hE / c
  double c, Z, time_e;
};

typedef enum {
  so_Riemann,
  so_Central,
  so_Alternating_sL,
  so_Alternating_sR,
  StarOption_SIZE
} StarOption;

typedef enum {
  cfem_1D,
  DG_2FUV,
  DG_1F_vStar,
  DG_1F_uStar,
  WeakFormulationType_SIZE
} WeakFormulationType;

typedef enum {
  dg_eps_m1 = -1,
  dg_eps_0,
  dg_eps_p1,
  DG_epsilonT_SIZE
} DG_epsilonT;

#endif
