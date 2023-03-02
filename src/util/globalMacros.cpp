#include "globalMacros.h"

// string g_prefileName = "";

Dcomplex Icomp(0.0, 1.0);
int serialNumber = 0;

#if DB_MODE
fstream db("debug.txt", ios::out);
#else
#define DB(x)
#endif

void setGlobalMembers()
{
    int x = 1;
}

