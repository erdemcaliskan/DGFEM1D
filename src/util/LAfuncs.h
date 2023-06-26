#ifndef LA_FUNCS__H
#define LA_FUNCS__H

#include "globalMacros.h"
#include <iostream>

class VECTOR;
class MATRIX;


typedef VECTOR DCVECTOR;
#if VCPP
typedef MATRIX DCMATRIX;
//typedef Eigen::MatrixXcd DCMATRIX;
//typedef Eigen::MatrixXf DCMATRIX;
//typedef Eigen::MatrixXd DCMATRIX;
#define ZEROMAT(x, n, m) { x.resize(n, m); x = 0.0;}
#define ZEROVEC(x, n) { x.resize(n); x = 0.0;}

#else // linux ...
#if USE_COMPLEX
#include <Eigen/Dense>
typedef Eigen::MatrixXcd DCMATRIX;
#define ZEROMAT(x, n, m)     x.setConstant(0.0)
#define ZEROVEC(x, n)     x.setConstant(0.0) //???
#else
#include <Eigen/Dense>
typedef Eigen::MatrixXd DCMATRIX;
#define ZEROMAT(x, n, m) x.setConstant(0.0)
#define ZEROVEC(x, n) { x.resize(n); x = 0.0;}
#endif
#endif

using namespace std;

#define USE_BLITZ 0

// Blitz is just one of the many libraries that provide vector and matrix
// storage capability. In C++ unlike Matlab there are no built-in matrix
// classes. Users include various matrix classes through including appropriate
// header files (*.h e.g. blitz/array.h) and perhaps calling provided libraries.
// These often do not come with standard C++ structures.

// In my computer I have installed blitz. To make it possible to do linear
// algebra matrix solution I have defined my own VECTOR and MATRIX classes
// below.
//  if USE_BLITZ is 1 blitz is used for matrix and vector. If not my definitions
//  (class VECTOR and MATRIX below) are used.

#if USE_BLITZ

#include "blitz/array.h"
typedef blitz::Array<double, 1> VECTOR;
typedef blitz::Array<double, 2> MATRIX;

#else

#include <vector> /* for vector */
using namespace std;

// C++ is a very powerful object-oriented language where as opposed to a process
// (subroutine) focuesed approach objects (instances of classes) interact.
// Classes refer to the set of data members and funcationalities that define a
// function (refer to the example below). For example Shape can be represented
// as a class with data members area and circumference and functions to read
// write the shape The example shape class will be discussed futher later on.

class VECTOR
{
    // the two functions below are operator overloads for output (<<) and input
    // (>>) operators in C++ built-in C++ members are printed and read by << and
    // >> User-defined classes need user defined functions so that << and >> would
    // make sense. User can specify what format an object is written or read as
    // can be seen in the definition of these function.

    // istream (input) and ostream (output) are c++ I/O streams (objects).

    // the "friend" keyword is often used with these function meaning that the
    // function can access "private" data members of the function. This makes
    // writing these functions easier as they often need to access all class
    // members.
    friend istream &operator>>(istream &in, VECTOR &dat);
    // const keywork means dat that is sent to the function cannot be changed.
    // It's a good programming practice to use the keyword const for data members
    // that are not supposed to be changed.
    friend ostream &operator<<(ostream &out, const VECTOR &dat);

    // keyword public means that functions and members after the keyword can be
    // access every where most class functions are public.
  public:
    VECTOR(unsigned int sizeIn = 0);
    void resize(unsigned int sizeIn);
    unsigned int size() const;
    unsigned int rows() const
    {
        return size();
    }
    // the following two are called "operator overloading"
    //		to enable direct use of paranthesis for indexing a matrix we
    // need to define operator()
    // then a(10) would make sense.

    // the first version can change the value and return the component of the
    // VECTOR at position i. A use of this is for example in a(10) = 15; (LHS)
    double &operator[](unsigned int i);
	double &operator()(unsigned int i);
	// second option is a "const" function. Adding const at the end of the
    // function name ensures that none of class members (vector<double> vec
    // herein) can be changed. Query functions are often suggested to be written
    // as const functions to prevent accidental changes to class data
    double operator[](unsigned int i) const;
	double operator()(unsigned int i) const;

    // This is another operator overloading that enables an = sign whose RHS is a
    // double. It enables setting the value for the entire VECTOR as we typically
    // write in mathematical notation. Otherwise something like
    //	VECTOR a(2);	// size 2
    // a = 10; would not make sense. Function below enables interpretation of the
    // expression in this line.
    VECTOR &operator=(double val);

    // functions and data after keyword "private"  are private and cannot be
    // access from outside of the class. generall class data and auxiliary
    // functions are private.
  private:
    // the following is a template. Type between < > says that this is a vector of
    // double. templates are commonly used in C++ as one class can represent
    // various types (e.g. vector of integers, doubles, ... by vecotr<int>,
    // vector<double>)
    vector<double> vec;
};

class MATRIX
{
    friend istream &operator>>(istream &in, MATRIX &dat);
    friend ostream &operator<<(ostream &out, const MATRIX &dat);

  public:
    MATRIX(unsigned int rowsIn = 0, unsigned int colsIn = 0);
    void resize(unsigned int rowsIn = 0, unsigned int colsIn = 0);
    void Multiply(const MATRIX &matIn, double factor = 1.0);
    unsigned int rows() const;
    unsigned int columns() const;
    vector<double> &operator[](unsigned int i);
    const vector<double> &operator[](unsigned int i) const;

	double& operator()(unsigned int i, unsigned int j);
	double operator()(unsigned int i, unsigned int j) const;
	double sum() const;

    //	double& operator()(unsigned int i, unsigned int j);
    //	double operator()(unsigned int i, unsigned int j) const;
    MATRIX &operator=(double val);

  private:
    vector<vector<double>> matx;
    int nrows;
    int ncols;
};
#endif

// this is a function that solves the system
// Ka = F. K and F are sent to the function and the result is written back to F
// (example in main.cpp) NOTE WE DO NOT INVERT K TO SOLVE FOR Ka = F. Rather we
// solve the system only for the given RHS F In the function below we use LU
// decomposition method. returns     singular (0 is good output 1 is not)
int LUsolve(MATRIX &K, VECTOR &F);

#endif