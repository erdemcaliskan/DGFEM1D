#ifndef LA_FUNCS_FINAL_STEP__H
#define LA_FUNCS_FINAL_STEP__H

#include "globalMacros.h"

void ReadV(VEC &vec, istream &in);
void PrintV(const VEC &vec, ostream &out);
void ReadM(MAT &mat, istream &in);
void PrintM(const MAT &mat, ostream &out);

void ReadVm1(VECm1 &vec, istream &in);
void PrintVm1(const VECm1 &vec, ostream &out);
void ReadMm1(MATm1 &mat, istream &in);
void PrintMm1(const MATm1 &mat, ostream &out);

void Inverse(const MAT &A, MAT &AInv);
#if DiM2a3_F
void Inverse(const MATm1 &A, MATm1 &AInv);
// condensedMat = mat_tt - mat_tn (mat_nn)^-1 mat_nt
void Break_nt_Matrix_n_t_components(bool computeInverse, const MAT &mat_nt,
                                    double &scalar_nn, VECm1 &vec_nt,
                                    VECm1 &vec_tn, MATm1 &mat_tt,
                                    MATm1 &mat_tt_condensed,
                                    MATm1 &mat_tt_condensed_Inv);
#endif

void CopyVec(const VECm1 &A, VECm1 &B);
void CopyVec(const VEC &A, VEC &B);
void CopyMat(const MATm1 &A, MATm1 &B);
void CopyMat(const MAT &A, MAT &B);
void CopyMat_withFactor(const MAT &A, MAT &B, double factor);
void setValue(VEC &A, double value = 0.0);
void setValue(MAT &A, double value = 0.0);
void setValue(VECm1 &A, double value = 0.0);
void setValue(MATm1 &A, double value = 0.0);
void ProductMatVec(const MAT &matA, const VEC &vecB, VEC &matA_times_vecB);
void ProductTransposeMatVec(const MAT &matA, const VEC &vecB,
                            VEC &matATranspose_times_vecB);
void ProductMatMat(const MAT &matA, const MAT &matB, MAT &matA_times_matB);
void AddVec(const VEC &vecA, const VEC &vecB, VEC &vecA_plus_vecB);
void SubtractVec(const VEC &vecA, const VEC &vecB, VEC &vecA_minus_vecB);
void FactorVec(VEC &vec, double factor);
void FactorMat(MAT &mat, double factor);
void LinearCombination2Vecs(const VEC &vecA, const VEC &vecB, double factorA,
                            double factorB, VEC &fAvecAplusfBvecB);

double Norm2(const VEC &vec);
double Normalize(const VEC &vec, VEC &normalized_vec);
// matrix is positive,
// row_i of eigenVecsRowWise and component i of eigenValues are eigenvector and
// eigenvalue #i since the matrix is positive, eigenvectors are orthonormal
// return is true if operation is successful
// eigenvalues >= 0 and are ordered from largest to smallest
bool EigenValuesEigenVectors_PositiveSymMatrix(const MAT &mat,
                                               MAT &eigenVecsRowWise,
                                               VEC &eigenValues);

#if DiM2a3_F
double Norm2m1(const Vc_dm1 &vec);
double Normalize_m1(const Vc_dm1 &vec, Vc_dm1 &normalized_vec);
void Extract_Shear(const VEC &vec_nt, Vc_dm1 &vec_t);
#endif
// vec = [v0, v1, v2] (3D), 2D and 1D smaller
// vec_n = v0, vec_n_pos = <vec_n>+
// Dim > 1	-> vec_t = [v1, v2], vec_t_abs = norm2(v1, v2), vec_effective =
// sqrt(vec_n_pos^2 + (beta * vec_t_abs)^2)
class Vec_nt {
public:
  Vec_nt();
  void SetZero();
  void Read(istream &in);
  void Compute_derivedValues_From_vec(const VEC &vecIn, double beta);
  inline double getEffective() const {
#if DiM1
    return vec_n_pos;
#else
    return vec_effective;
#endif
  }
  // input
  VEC vec;
  //

  // computed
  bool nt_processed;
  double vec_n;
  double vec_n_pos;
#if DiM2a3_F
  Vc_dm1 vec_t;
  double vec_t_abs;

private:
  double vec_effective;
#endif
private:
  void Compute_derivedValues_From_vec(double beta);
};

template <class T> void ReadVector(istream &in, vector<T> &dat) {
  T tmp;
  string buf;
  // Error: Returning pointer to local variable 'tmpstr' that will be invalid
  // when returning.
  READ_NSTRING(in, buf, buf);
  dat.clear();
  if (buf != "{")
    THROW("istream should start with {");
  READ_NSTRING(in, buf, buf);
  while (buf != "}") {
    if (fromString(buf, tmp) == false) {
      cout << "buf\t" << buf << '\n';
      THROW("Value must be integer\n");
    }
    dat.push_back(tmp);
    READ_NSTRING(in, buf, buf);
  }
}

#endif
