#ifndef GLOBAL_MACROS__H
#define GLOBAL_MACROS__H

#include "LAfuncsFixed.h"
#include "commonMacros.h"
#include <complex>

typedef Vc_dd VEC;
typedef Mtrx_dd MAT;

#if 1 // DiM2a3_F
typedef Vc_dm1 VECm1;
typedef Mtrx_dm1 MATm1;
#endif

#ifndef EXIT
#define EXIT                                                                   \
  exit(1);                                                                     \
  getchar();                                                                   \
  getchar();
#endif
#ifndef THROW
#define THROW(msg)                                                             \
  {                                                                            \
    {                                                                          \
      char tmpstr[255];                                                        \
      sprintf(tmpstr, "In %s, line %d : %s \n", __FILE__, __LINE__, msg);      \
      cerr << tmpstr;                                                          \
      getchar();                                                               \
      getchar();                                                               \
      throw(tmpstr);                                                           \
    }                                                                          \
  }
#endif

#define Dcomplex std::complex<double>
#define NUMBR Dcomplex

// constexpr auto NUM_SIDES = 2;
// constexpr auto SDL = 0;
// constexpr auto SDR = 1;

typedef long GID;

// extern string g_prefileName;
extern Dcomplex Icomp;

void setGlobalMembers();

#endif