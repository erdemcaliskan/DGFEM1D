#ifndef DIMS__H
#define DIMS__H

#include <array>
#include <assert.h> /* assert */
#include <cstdlib>
#include <cstring>
#include <fstream>  /* file streams for input and output */
#include <iostream> /* for input and output using >> and << operator */
#include <math.h>   /* math functions */
#include <sstream>
#include <string>

#define USE_DEALII 0
#define USE_DEALII_VECMATIO 0	//1

using namespace std;

#define DiM1 1
#define DiM2 0
#define DiM3 0

#if DiM2
#define DiM 2
#define DiMm1 1
#define DiM2a3(x) x
#define DiM2a3_F 1
#else
#if DiM3
#define DiM 3
#define DiMm1 2
#define DiM2a3(x) x
#define DiM2a3_F 1
#else
#define DiM 1
#define DiMm1 0
#define DiM2a3(x)
#define DiM2a3_F 0
#endif
#endif

#define MAX(a, b) ((a) > (b) ? a : b)
#define MIN(a, b) ((a) > (b) ? b : a)

#ifndef EXIT
#define EXIT                                                                                                           \
    exit(1);                                                                                                           \
    getchar();                                                                                                         \
    getchar();
#endif

// to exit the code with useful information when something is wrong. It provides
// file and line number with a message if a check is incorrect. example
//					if (denominator == 0) THROW("cannot
// divide be zero!\n");
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

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

constexpr auto NUM_SIDES = 2;
constexpr auto SDL = 0;
constexpr auto SDR = 1;

#endif
