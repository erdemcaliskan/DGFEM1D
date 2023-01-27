#ifndef ONE_DIMENSIONAL_ELEMENT__H
#define ONE_DIMENSIONAL_ELEMENT__H

#include "LAFuncs.h"
#include "globalTypesClasses.h"

// this is the parent element from xi = -1 to 1 and no real geometry or material properties are considered here
// N is shape function
// Bxi = dN / dxi
class OneDimensionalParentElement
{
  public:
    // element order
    unsigned int polyOrder;
    // whether we want the mass-lumpting option or not
    bool lumpMass;

    void Initialize(unsigned int polyOrderIn, bool lumpMassIn);
    void CalculateN(double xi, VECTOR &N) const; // N resized ndof and stores shape functions
    void CalculateB_xi(double xi,
                       VECTOR &Bxi) const; // Bxi resized ndof and stores derivative shape functions w.r.t. xi
    /// computed
    unsigned int ndof;
    MATRIX mpe, kpe; // me = int_xi=-1 -> 1 N^t N d xi,   ke = int_xi=-1 -> 1 Bxi^t Bxi d xi

    VECTOR Ne_leftNode, Ne_rightNode, Be_leftNode, Be_rightNode;
};

class OneDimensionalElement
{
  public:
    double GetJacobian() const; // this is a constant Jacobian element (element is not distorted)
    double x_2xi(double x, double xLeft) const;
    void CalculateB(const OneDimensionalParentElement &parent, double xi,
                    VECTOR &B) const; // B resized ndof and stores derivative shape functions w.r.t. x (real coordinate)
    void CalculateB(const OneDimensionalParentElement &parent, double x, double xLeft, VECTOR &B) const;
    // CalculateN from (xi) directly uses CalculateN from the parent element
    // CalculateN from (x) similar to above x -> xi, then call parent calculate N

    ElementProperties elementProps;
    MATRIX me, ce, ke;
};

#endif