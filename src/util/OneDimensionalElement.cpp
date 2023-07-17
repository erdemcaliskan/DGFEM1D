#include "OneDimensionalElement.h"

void OneDimensionalParentElement::Initialize(unsigned int polyOrderIn, bool lumpMassIn)
{
    lumpMass = lumpMassIn;
    polyOrder = polyOrderIn;
    ndof = polyOrder + 1;
    kpe.resize(ndof, ndof);
    mpe.resize(ndof, ndof);

    Ne_leftNode.resize(ndof);
    Ne_rightNode.resize(ndof);
    Be_leftNode.resize(ndof);
    Be_rightNode.resize(ndof);
    Ne_leftNode = 0.0;
    Ne_rightNode = 0.0;
    Ne_leftNode[0] = 1; // shape function 0 is 1 at the left node
    if (polyOrder > 0)
        Ne_rightNode[1] = 1; // shape function 0 is 1 at the right node
    Be_leftNode = 0.0;
    Be_rightNode = 0.0;

    if (polyOrder == 0)
    {
        kpe = 0.0;
        mpe[0][0] = 2.0;
        Ne_leftNode[0] = 1;
        Ne_rightNode[0] = 1;
    }
    else if (polyOrder == 1)
    {
        if (!lumpMass)
        {
            double fact = 1.0 / 3.0;
            mpe[0][0] = mpe[1][1] = 2.0 * fact;
            mpe[0][1] = fact;
            mpe[1][0] = fact;
        }
        else
        {
            mpe[0][0] = mpe[1][1] = 1.0;
            mpe[0][1] = 0.0;
            mpe[1][0] = 0.0;
        }
        // stiffness
        kpe[0][0] = kpe[1][1] = 0.5;
        kpe[0][1] = -0.5;
        kpe[1][0] = -0.5;

        // dN/dxi on the left side
        Be_leftNode[0] = -0.5;
        Be_leftNode[1] = 0.5;

        // dN/dxi on the right side
        Be_rightNode[0] = -0.5;
        Be_rightNode[1] = 0.5;
    }
    else if (polyOrder == 2)
    {
        // check these mass matrices and stiffness
        /// EDITED
        if (!lumpMass)
        {
            double fact = 1.0 / 15.0;
            mpe[0][0] = mpe[1][1] = 4.0 * fact;
            mpe[0][1] = mpe[1][0] = -1.0 * fact;
            mpe[2][2] = 16.0 * fact;
            mpe[0][2] = mpe[2][0] = mpe[1][2] = mpe[2][1] = 2.0 * fact;
        }
        else
        {
            double fact = 1.0 / 3.0;
            mpe[0][0] = mpe[1][1] = fact;
            mpe[0][1] = mpe[1][0] = 0.0;
            mpe[2][2] = 4.0 * fact; // 2023/06/24 was 2.0 * fact;
            mpe[0][2] = mpe[2][0] = mpe[1][2] = mpe[2][1] = 0.0;
        }
        // stiffness
        /// EDITED
        double fact = 1.0 / 6.0;
        kpe[0][0] = kpe[1][1] = 7.0 * fact;
        kpe[0][1] = kpe[1][0] = fact;
        kpe[2][2] = 16.0 * fact;
        kpe[0][2] = kpe[2][0] = kpe[1][2] = kpe[2][1] = -8.0 * fact;

        // dN/dxi on the left side
        Be_leftNode[0] = -1.5;
        Be_leftNode[1] = -0.5;
        Be_leftNode[2] = 2.0;

        // dN/dxi on the right side
        Be_rightNode[0] = 0.5;
        Be_rightNode[1] = 1.5;
        Be_rightNode[2] = -2.0;
    }
    else if (polyOrder == 3)
    {
        // check these mass matrices and stiffness
        /// EDITED
        if (!lumpMass)
        {
            double fact = 1.0 / 68040.0;
            mpe[0][0] = mpe[1][1] = 128.0 * fact;
            mpe[0][1] = mpe[1][0] = 152.0 * fact;
            mpe[0][2] = mpe[2][0] = -12.0 * fact;
            mpe[0][3] = mpe[3][0] = 33.0 * fact;
            mpe[2][2] = mpe[3][3] = 72.0 * fact;
            mpe[1][2] = mpe[2][1] = 33.0 * fact;
            mpe[1][3] = mpe[3][1] = -12.0 * fact;
            mpe[2][3] = mpe[3][2] = -9.0 * fact;
        }
        else
        {
            double fact = 1.0 / 9720.0;
            mpe[0][0] = mpe[1][1] = 43.0 * fact;
            mpe[0][1] = mpe[1][0] = 0.0;
            mpe[0][2] = mpe[2][0] = 0.0;
            mpe[0][3] = mpe[3][0] = 0.0;
            mpe[2][2] = mpe[3][3] = 12.0 * fact;
            mpe[1][2] = mpe[2][1] = 0.0;
            mpe[1][3] = mpe[3][1] = 0.0;
            mpe[2][3] = mpe[3][2] = 0.0;
        }
        // stiffness
        /// EDITED
        double fact = 1.0 / 6480.0;
        kpe[0][0] = kpe[1][1] = 148.0 * fact;
        kpe[0][1] = kpe[1][0] = -13.0 * fact;
        kpe[0][2] = kpe[2][0] = 18.0 * fact;
        kpe[0][3] = kpe[3][0] = -63.0 * fact;
        kpe[2][2] = kpe[3][3] = 48.0 * fact;
        kpe[1][2] = kpe[2][1] = -63.0 * fact;
        kpe[1][3] = kpe[3][1] = 18.0 * fact;
        kpe[2][3] = kpe[3][2] = -33.0 * fact;

        // dN/dxi on the left side
        Be_leftNode[0] = -11.0 / 36.0;
        Be_leftNode[1] = 1.0 / 18.0;
        Be_leftNode[2] = -1.0 / 12.0;
        Be_leftNode[3] = 1.0 / 6.0;  

        // dN/dxi on the right side
        Be_rightNode[0] = -1.0 / 18.0;
        Be_rightNode[1] = 11.0 / 36.0;
        Be_rightNode[2] = -1.0 / 6.0;  
        Be_rightNode[3] = 1.0 / 12.0;
    }
}

void OneDimensionalParentElement::CalculateN(double xi, VECTOR &N) const
{
    THROW("complete the function\n");
}

void OneDimensionalParentElement::CalculateB_xi(double xi, VECTOR &Bxi) const
{
    THROW("complete the function\n");
}

double OneDimensionalElement::GetJacobian() const
{
    return 0.5 * elementProps.hE;
}

void OneDimensionalElement::CalculateB(const OneDimensionalParentElement &parent, double xi, VECTOR &B) const
{
    parent.CalculateB_xi(xi, B);
    double Jacobian = GetJacobian();
    for (unsigned int i = 0; i < parent.ndof; ++i)
        B[i] *= Jacobian;
}

double OneDimensionalElement::x_2xi(double x, double xLeft) const
{
    // x = N0Linear(xi) xL + N1Linear(xi) XR = (1 - xi)/2 xL + (1 + xi/2) xR = (1 - xi)/2 (xL) + (1 + xi/s)(xL + he) =
    // xL + (1 + xi)/2 * he
    double xi = 2.0 * (x - xLeft) / elementProps.hE - 1.0;
    return xi;
}

void OneDimensionalElement::CalculateB(const OneDimensionalParentElement &parent, double x, double xLeft,
                                       VECTOR &B) const
{
    double xi = x_2xi(x, xLeft);
    CalculateB(parent, xi, B);
}