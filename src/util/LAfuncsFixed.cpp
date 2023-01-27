#include "LAfuncsFixed.h"

ostream &operator<<(ostream &out, const Vc_dd &dat)
{
#if DiM1
    out << dat[0];
#endif
#if DiM2
    out << dat[0] << '\t' << dat[1];
#endif
#if DiM3
    out << dat[0] << '\t' << dat[1] << '\t' << dat[2];
#endif
    return out;
}

istream &operator>>(istream &in, Vc_dd &dat)
{
#if DiM1
    in >> dat[0];
#endif
#if DiM2
    in >> dat[0] >> dat[1];
#endif
#if DiM3
    in >> dat[0] >> dat[1] >> dat[2];
#endif
    return in;
}

ostream &operator<<(ostream &out, const Mtrx_dd &dat)
{
#if DiM1
    out << dat[0][0];
#endif
#if DiM2
    out << dat[0][0] << '\t' << dat[0][1] << '\n';
    out << dat[1][0] << '\t' << dat[1][1];
#endif
#if DiM3
    out << dat[0][0] << '\t' << dat[0][1] << '\t' << dat[0][2] << '\n';
    out << dat[1][0] << '\t' << dat[1][1] << '\t' << dat[1][2] << '\n';
    out << dat[2][0] << '\t' << dat[2][1] << '\t' << dat[2][2];
#endif
    return out;
}

istream &operator>>(istream &in, Mtrx_dd &dat)
{
#if DiM1
    in >> dat[0][0];
#endif
#if DiM2
    in >> dat[0][0] >> dat[0][1];
    in >> dat[1][0] >> dat[1][1];
#endif
#if DiM3
    in >> dat[0][0] >> dat[0][1] >> dat[0][2];
    in >> dat[1][0] >> dat[1][1] >> dat[1][2];
    in >> dat[2][0] >> dat[2][1] >> dat[2][2];
#endif
    return in;
}

#if 1 // DiM2a3_F

ostream &operator<<(ostream &out, const Vc_dm1 &dat)
{
#if DiM2
    out << dat[0];
#endif
#if DiM3
    out << dat[0] << '\t' << dat[1];
#endif
    return out;
}

istream &operator>>(istream &in, Vc_dm1 &dat)
{
#if DiM2
    in >> dat[0];
#endif
#if DiM3
    in >> dat[0] >> dat[1];
#endif
    return in;
}

ostream &operator<<(ostream &out, const Mtrx_dm1 &dat)
{
#if DiM2
    out << dat[0][0];
#endif
#if DiM3
    out << dat[0][0] << '\t' << dat[0][1] << '\n';
    out << dat[1][0] << '\t' << dat[1][1];
#endif
    return out;
}

istream &operator>>(istream &in, Mtrx_dm1 &dat)
{
#if DiM2
    in >> dat[0][0];
#endif
#if DiM3
    in >> dat[0][0] >> dat[0][1];
    in >> dat[1][0] >> dat[1][1];
#endif
    return in;
}

#endif
