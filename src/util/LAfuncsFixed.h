#ifndef LA_FUNCS_FIXED__H
#define LA_FUNCS_FIXED__H

#include <iostream>
using namespace std;
#include "Dims.h"

typedef double Vc_dd[DiM];
typedef double Mtrx_dd[DiM][DiM];

ostream &operator<<(ostream &out, const Vc_dd &dat);
istream &operator>>(istream &in, Vc_dd &dat);
ostream &operator<<(ostream &out, const Mtrx_dd &dat);
istream &operator>>(istream &in, Mtrx_dd &dat);

#if DiM2a3_F
typedef double Vc_dm1[DiMm1];
typedef double Mtrx_dm1[DiMm1][DiMm1];
#else
typedef double Vc_dm1[2];
typedef double Mtrx_dm1[2][2];
#endif

ostream &operator<<(ostream &out, const Vc_dm1 &dat);
istream &operator>>(istream &in, Vc_dm1 &dat);
ostream &operator<<(ostream &out, const Mtrx_dm1 &dat);
istream &operator>>(istream &in, Mtrx_dm1 &dat);

#endif
