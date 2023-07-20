#ifndef LA_FUNCS__H
#define LA_FUNCS__H

#include "globalMacros.h"
#include <iostream>
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

#if	USE_BLITZ
#else

#include <vector> /* for vector */
using namespace std;

// C++ is a very powerful object-oriented language where as opposed to a process
// (subroutine) focuesed approach objects (instances of classes) interact.
// Classes refer to the set of data members and funcationalities that define a
// function (refer to the example below). For example Shape can be represented
// as a class with data members area and circumference and functions to read
// write the shape The example shape class will be discussed futher later on.

template <class T>
class rVEC
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
	friend istream &operator>>(istream &in, rVEC<T> &dat)
	{
		string buf;
		unsigned int sz;
#if USE_DEALII_VECMATIO
		sz = DiM;
#else
		in >> buf >> sz;
#endif
		dat.resize(sz);
		for (unsigned int i = 0; i < sz; ++i)
			in >> dat[i];
		return in;
	}
    // const keywork means dat that is sent to the function cannot be changed.
    // It's a good programming practice to use the keyword const for data members
    // that are not supposed to be changed.
	friend ostream &operator<<(ostream &out, const rVEC<T> &dat)
	{
		unsigned int sz = dat.size();
#if !USE_DEALII_VECMATIO
		out << "size\t" << sz << '\t';
#endif
		for (unsigned int i = 0; i < sz; ++i)
			out << dat[i] << '\t';
		return out;
	}

    // keyword public means that functions and members after the keyword can be
    // access every where most class functions are public.
  public:
    rVEC<T>(unsigned int sizeIn = 0);
	void resize(unsigned int sizeIn);
	unsigned int size() const;
    inline unsigned int rows() const    {        return size();    }
    // the following two are called "operator overloading"
    //		to enable direct use of paranthesis for indexing a matrix we
    // need to define operator()
    // then a(10) would make sense.

    // the first version can change the value and return the component of the
    // VECTOR at position i. A use of this is for example in a(10) = 15; (LHS)
	inline T &operator[](unsigned int i) { return vec[i]; }
	inline T &operator()(unsigned int i) { return vec[i]; }
	// second option is a "const" function. Adding const at the end of the
    // herein) can be changed. Query functions are often suggested to be written
    // as const functions to prevent accidental changes to class data
	inline T operator[](unsigned int i) const { return vec[i]; }
	inline T operator()(unsigned int i) const { return vec[i]; }

    // This is another operator overloading that enables an = sign whose RHS is a
    // . It enables setting the value for the entire VECTOR as we typically
    // write in mathematical notation. Otherwise something like
    //	VECTOR a(2);	// size 2
    // a = 10; would not make sense. Function below enables interpretation of the
    // expression in this line.
	rVEC<T> &operator=(T val);
	rVEC<T> &operator=(rVEC<T>& other);

    // functions and data after keyword "private"  are private and cannot be
    // access from outside of the class. generall class data and auxiliary
    // functions are private.
  private:
    // the following is a template. Type between < > says that this is a vector of
    // T. templates are commonly used in C++ as one class can represent
    // various types (e.g. vector of integers, Ts, ... by vecotr<int>,
    // vector<T>)
    vector<T> vec;
};

template <class T>
class rMAT
{
    friend istream &operator>>(istream &in, rMAT<T> &dat)
	{
		string buf;
		unsigned int nrow, ncol;
#if USE_DEALII_VECMATIO
		nrow = DiM, ncol = DiM;
#else
		in >> buf >> nrow >> buf >> ncol;
#endif
		dat.resize(nrow, ncol);
		for (unsigned int i = 0; i < nrow; ++i)
		{
			for (unsigned int j = 0; j < ncol; ++j)
				in >> dat[i][j];
		}
		return in;
	}

    friend ostream &operator<<(ostream &out, const rMAT<T> &dat)
	{
		unsigned int nrow = dat.rows(), ncol = dat.columns();
#if !USE_DEALII_VECMATIO
		out << "rows\t" << nrow << "\tcols\t" << ncol << '\n';
#endif
		for (unsigned int i = 0; i < nrow; ++i)
		{
			for (unsigned int j = 0; j < ncol; ++j)
				out << dat[i][j] << '\t';
#if USE_DEALII_VECMATIO
			out << '\t';
#else
			out << '\n';
#endif
		}
		return out;
	}

  public:
    rMAT<T>(unsigned int rowsIn = 0, unsigned int colsIn = 0);
    void resize(unsigned int rowsIn = 0, unsigned int colsIn = 0);
    void Multiply(const rMAT<T> &matIn, T factor);
	inline unsigned int rows() const { return nrows;	};
    unsigned int columns() const {	return ncols;	}

	vector<T> &operator[](unsigned int i) { return matx[i]; };
    const vector<T> &operator[](unsigned int i) const {		return matx[i];	}

	T& operator()(unsigned int i, unsigned int j)	{		return matx[i][j];	}

	T operator()(unsigned int i, unsigned int j) const { return matx[i][j]; }
	T sum() const;

    rMAT<T> &operator=(T val);
	rMAT<T> &operator=(rMAT<T>& matOther);

  private:
    vector <vector<T> > matx;
    int nrows;
    int ncols;
};
#endif

// this is a function that solves the system
// Ka = F. K and F are sent to the function and the result is written back to F
// (example in main.cpp) NOTE WE DO NOT INVERT K TO SOLVE FOR Ka = F. Rather we
// solve the system only for the given RHS F In the function below we use LU
// decomposition method. returns     singular (0 is good output 1 is not)
template <class T>
int LUsolve(rMAT<T> &K, rVEC<T> &F);


template <class T>
rVEC<T>::rVEC(unsigned int sizeIn)
{
	resize(sizeIn);
}

template <class T>
void rVEC<T>::resize(unsigned int sizeIn)
{
	vec.resize(sizeIn);
}

template <class T>
unsigned int rVEC<T>::size() const
{
	return vec.size();
}

template <class T>
rVEC<T>& rVEC<T>::operator=(T val)
{
	unsigned int sz = vec.size();
	for (unsigned int i = 0; i < sz; ++i)
		vec[i] = val;
	return *this;
}

template <class T>
rVEC<T>& rVEC<T>::operator=(rVEC<T>& other)
{
	unsigned int sz = other.vec.size();
	vec.resize(sz);
	for (unsigned int i = 0; i < sz; ++i)
		vec[i] = other.vec[i];
	return *this;
}

template <class T>
rMAT<T>::rMAT(unsigned int rowsIn, unsigned int colsIn)
{
	resize(rowsIn, colsIn);
}

template <class T>
void rMAT<T>::resize(unsigned int rowsIn, unsigned int colsIn)
{
	nrows = rowsIn;
	ncols = colsIn;
	if (nrows < 0)
	{
		cout << "nrows\t" << nrows << '\n';
		THROW("nrows < 0");
	}
	matx.resize(rowsIn);
	for (unsigned int i = 0; i < nrows; ++i)
	{
		matx[i].resize(ncols);
	}
}

template <class T>
void rMAT<T>::Multiply(const rMAT<T> &matIn, T factor)
{
	nrows = matIn.nrows;
	ncols = matIn.ncols;

	matx.resize(nrows);
	for (unsigned int i = 0; i < nrows; ++i)
	{
		matx[i].resize(ncols);
		for (unsigned int j = 0; j < ncols; ++j)
			matx[i][j] = factor * matIn[i][j];
	}
}

template <class T>
T rMAT<T>::sum() const
{
	T sm = 0;
	for (unsigned int i = 0; i < nrows; ++i)
		for (unsigned int j = 0; j < ncols; ++j)
			sm += matx[i][j];
	return sm;
}

template <class T>
rMAT<T>& rMAT<T>::operator=(T val)
{
	for (unsigned int i = 0; i < nrows; ++i)
	{
		for (unsigned int j = 0; j < ncols; ++j)
			matx[i][j] = val;
	}
	return *this;
}

template <class T>
rMAT<T>& rMAT<T>::operator=(rMAT<T>& matOther)
{
	resize(matOther.nrows, matOther.ncols);
	for (unsigned int i = 0; i < nrows; ++i)
	{
		for (unsigned int j = 0; j < ncols; ++j)
			matx[i][j] = matOther.matx[i][j];
	}
	return *this;
}

#if USE_BLITZ
#include "blitz/array.h"
typedef blitz::array<double, 1> VECTOR;
typedef blitz::array<Dcomplex, 1> CVECTOR;
typedef blitz::array<double, 2> MATRIX;
#else
typedef rVEC<double> VECTOR;
typedef rVEC<Dcomplex> CVECTOR;
typedef rMAT<double> MATRIX;
#endif

#if VCPP
#if USE_BLITZ
typedef blitz::array<Dcomplex, 2> CMATRIX;
#else
typedef rMAT<Dcomplex> CMATRIX;
#endif

#if USE_COMPLEX
typedef CVECTOR DCVECTOR;
typedef CMATRIX DCMATRIX;
#else
typedef VECTOR DCVECTOR;
typedef MATRIX DCMATRIX;
#endif

#define ZEROMAT(x, n, m) { x.resize(n, m); x = 0;}
#define ZEROVEC(x, n) { x.resize(n); x = 0;}

#else // linux ...
#if USE_COMPLEX
typedef CVECTOR DCVECTOR;
#include <Eigen/Dense>
typedef Eigen::MatrixXcd DCMATRIX;
#define ZEROMAT(x, n, m)     x.setConstant(0)
#define ZEROVEC(x, n)     x.setConstant(0) //???
#else
#include <Eigen/Dense>
typedef VECTOR DCVECTOR;
typedef Eigen::MatrixXd DCMATRIX;
#define ZEROMAT(x, n, m) x.setConstant(0)
#define ZEROVEC(x, n) { x.resize(n); x = 0;}
#endif
#endif

#endif

// 	This is the LU decomposition and solve implementation that we used
// 	with TNT.  I's put here so people don't need to use CLAPACK to run
// 	the code.
template <class T>
int LUsolve(rMAT<T> &K, rVEC<T> &F)
{
	int m, n, pivsign;
	m = K.rows();
	n = K.columns();
	rMAT<T> LU_(m, n);
	LU_ = K; // CHECK: can we just get along with referencing? It seems to be the
			 // case, then we save some mempory copy here
			 //    MATRIX LU_(K); // refernce option

	int i = 0;
	int j = 0;
	int k = 0;

	pivsign = 1;
	rVEC<T> LUcolj(m);
	vector<int> piv(m);

	for (i = 0; i < m; i++)
		piv[i] = i;

	for (j = 0; j < n; j++)
	{
		// Make a copy of the j-th column to localize references.
		for (i = 0; i < m; i++)
			LUcolj[i] = LU_[i][j];

		// Apply previous transformations.
		for (i = 0; i < m; i++)
		{
			int kmax = i < j ? i : j; // min(i,j);
			T s = 0;
			for (k = 0; k < kmax; k++)
				s += LU_[i][k] * LUcolj[k];

			LU_[i][j] = LUcolj[i] -= s;
		}

		// Find pivot and exchange if necessary.
		int p = j;
		for (int i = j + 1; i < m; i++)
#if USE_COMPLEX
			if (abs(LUcolj[i]) > abs(LUcolj[p]))
#else
			if (fabs(LUcolj[i]) > fabs(LUcolj[p]))
#endif
				p = i;

		if (p != j)
		{
			for (k = 0; k < n; k++)
			{
				T t = LU_[p][k];
				LU_[p][k] = LU_[j][k];
				LU_[j][k] = t;
			}
			k = piv[p];
			piv[p] = piv[j];
			piv[j] = k;
			pivsign = -pivsign;
		}

		// Compute multipliers.
		if ((j < m) && (LU_[j][j] != 0.0))
			for (int i = j + 1; i < m; i++) //... [thite, 6/30/2003]
				LU_[i][j] /= LU_[j][j];
	}

	int issingular = 0;
	for (int j = 0; j < n; j++)
	{
		if (LU_[j][j] == 0.0)
			issingular = 1;
	}
#if 0
	if (issingular != 0)
	{
		std::cout << "Solve procedure failed:  FILE:  " << __FILE__

			<< ", LINE:  " << __LINE__ << std::endl;
		return issingular;
		//        exit(1);
	}
#endif
	int piv_length = piv.size();
	DCVECTOR x(piv_length);
	for (int I = 0; I < piv_length; I++)
		x[I] = F[piv[I]];

	// Solve L*Y = B(piv)
	for (k = 0; k < n; k++)
		for (i = k + 1; i < n; i++)
			x[i] -= x[k] * LU_[i][k];

	// Solve U*X = Y;
	for (k = n - 1; k >= 0; k--)
	{
		x[k] /= LU_[k][k];
		for (i = 0; i < k; i++)
			x[i] -= x[k] * LU_[i][k];
	}

	F = x;
	K = LU_;
	return issingular;
} // end TNT version of LUsolve with 1D F