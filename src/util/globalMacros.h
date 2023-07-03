#ifndef GLOBAL_MACROS__H
#define GLOBAL_MACROS__H

#include "LAfuncsFixed.h"
#include "commonMacros.h"
#include <complex>

#if VCPP
#define USE_COMPLEX	0
//#define USE_COMPLEX	1
#else
//#define USE_COMPLEX	1
#define USE_COMPLEX	0
#endif

typedef Vc_dd VEC;
typedef Mtrx_dd MAT;

#if 1 // DiM2a3_F
typedef Vc_dm1 VECm1;
typedef Mtrx_dm1 MATm1;
#endif

#ifndef EXIT
#define EXIT                                                                                                           \
    exit(1);                                                                                                           \
    getchar();                                                                                                         \
    getchar();
#endif
#ifndef THROW
#define THROW(msg)                                                                                                     \
    {                                                                                                                  \
        {                                                                                                              \
            char tmpstr[255];                                                                                          \
            sprintf(tmpstr, "In %s, line %d : %s \n", __FILE__, __LINE__, msg);                                        \
            cerr << tmpstr;                                                                                            \
            getchar();                                                                                                 \
            getchar();                                                                                                 \
            throw(tmpstr);                                                                                             \
        }                                                                                                              \
    }
#endif

#define Dcomplex std::complex<double>
#if USE_COMPLEX
#define NUMBR Dcomplex
#else
#define NUMBR double
#endif

// constexpr auto NUM_SIDES = 2;
// constexpr auto SDL = 0;
// constexpr auto SDR = 1;

typedef long GID;

// extern string g_prefileName;
extern Dcomplex Icomp;
extern int serialNumber;
void setGlobalMembers();

#define DB_MODE 1

#if DB_MODE
#define DB(x) x
extern fstream db;
#else
#define DB(x)
#endif


#endif