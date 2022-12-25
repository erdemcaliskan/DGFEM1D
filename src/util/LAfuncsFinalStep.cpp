#include "LAfuncsFinalStep.h"

void ReadV(VEC &vec, istream &in) {
  string buf;
  in >> vec[0];
#if DiM2a3_F
  in >> vec[1];
#if DiM3
  in >> vec[2];
#endif
#endif
}

void PrintV(const VEC &vec, ostream &out) {
  out << vec[0];
#if DiM2a3_F
  out << '\t' << vec[1];
#if DiM3
  out << '\t' << vec[2];
#endif
#endif
}

void ReadM(MAT &mat, istream &in) {
#if DiM1
  in >> mat[0][0];
#endif
#if DiM2
  in >> mat[0][0] >> mat[0][1];
  in >> mat[1][0] >> mat[1][1];
#endif
#if DiM3
  in >> mat[0][0] >> mat[0][1] >> mat[0][2];
  in >> mat[1][0] >> mat[1][1] >> mat[1][2];
  in >> mat[2][0] >> mat[2][1] >> mat[2][2];
#endif
}

void PrintM(const MAT &mat, ostream &out) {
#if DiM1
  out << mat[0][0];
#endif
#if DiM2
  out << '\t' << mat[0][0] << '\t' << mat[0][1] << '\n';
  out << '\t' << mat[1][0] << '\t' << mat[1][1];
#endif
#if DiM3
  out << '\t' << mat[0][0] << '\t' << mat[0][1] << '\t' << mat[0][2] << '\n';
  out << '\t' << mat[1][0] << '\t' << mat[1][1] << '\t' << mat[1][2] << '\n';
  out << '\t' << mat[2][0] << '\t' << mat[2][1] << '\t' << mat[2][2];
#endif
}

void ReadVm1(VECm1 &vec, istream &in) {
#if DiM2a3_F
  in >> vec[0];
#if DiM3
  in >> vec[1];
#endif
#endif
}

void PrintVm1(const VECm1 &vec, ostream &out) {
#if DiM2a3_F
  out << '\t' << vec[0];
#if DiM3
  out << '\t' << vec[1];
#endif
#endif
}

void ReadMm1(MATm1 &mat, istream &in) {
#if DiM2
  in >> mat[0][0];
#endif
#if DiM3
  in >> mat[0][0] >> mat[0][1];
  in >> mat[1][0] >> mat[1][1];
#endif
}

void PrintMm1(const MATm1 &mat, ostream &out) {
#if DiM2
  out << mat[0][0];
#endif
#if DiM3
  out << '\t' << mat[0][0] << '\t' << mat[0][1] << '\n';
  out << '\t' << mat[1][0] << '\t' << mat[1][1];
#endif
}

void Inverse(const MAT &A, MAT &AInv) {
#if DiM1
  AInv[0][0] = 1.0 / A[0][0];
  return;
#endif
#if DiM2
  double det = A[1][1] * A[0][0] - A[0][1] * A[1][0];
  double detInv = 1.0 / det;
  AInv[0][0] = detInv * A[1][1];
  AInv[1][1] = detInv * A[0][0];
  AInv[0][1] = -detInv * A[0][1];
  AInv[1][0] = -detInv * A[1][0];
  return;
#endif
  THROW("Inverse not implemented for D3\n");
}

#if DiM2a3_F
void Inverse(const MATm1 &A, MATm1 &AInv) {
#if DiM2
  AInv[0][0] = 1.0 / A[0][0];
  return;
#endif
#if DiM3
  double det = A[1][1] * A[0][0] - A[0][1] * A[1][0];
  double detInv = 1.0 / det;
  AInv[0][0] = detInv * A[1][1];
  AInv[1][1] = detInv * A[0][0];
  AInv[0][1] = -detInv * A[0][1];
  AInv[1][0] = -detInv * A[1][0];
  return;
#endif
}

void Break_nt_Matrix_n_t_components(bool computeInverse, const MAT &mat_nt,
                                    double &scalar_nn, VECm1 &vec_nt,
                                    VECm1 &vec_tn, MATm1 &mat_tt,
                                    MATm1 &mat_tt_condensed,
                                    MATm1 &mat_tt_condensed_Inv) {
  scalar_nn = mat_nt[0][0];
  int ip1, jp1;
  double tmp;
  for (int i = 0; i < DiMm1; ++i) {
    ip1 = i + 1;
    vec_nt[i] = mat_nt[0][ip1];
    vec_tn[i] = mat_nt[ip1][0];
    for (int j = 0; j < DiMm1; ++j) {
      jp1 = j + 1;
      tmp = mat_nt[ip1][jp1];
      mat_tt[i][j] = tmp;
      mat_tt_condensed[i][j] = tmp;
    }
  }
  tmp = 1.0 / scalar_nn;
  for (int i = 0; i < DiMm1; ++i) {
    for (int j = 0; j < DiMm1; ++j)
      mat_tt_condensed[i][j] -= vec_tn[i] * tmp * vec_nt[j];
  }
  if (computeInverse)
    Inverse(mat_tt_condensed, mat_tt_condensed_Inv);
}
#endif

void CopyVec(const VECm1 &A, VECm1 &B) {
#if DiM2
  B[0] = A[0];
#endif
#if DiM3
  B[0] = A[0];
  B[1] = A[1];
#endif
}

void CopyVec(const VEC &A, VEC &B) {
  for (int i = 0; i < DiM; ++i)
    B[i] = A[i];
}

void CopyMat(const MATm1 &A, MATm1 &B) {
#if DiM2
  B[0][0] = A[0][0];
#endif
#if DiM3
  B[0][0] = A[0][0];
  B[0][1] = A[0][1];
  B[1][0] = A[1][0];
  B[1][1] = A[1][1];
#endif
}

void CopyMat(const MAT &A, MAT &B) {
  for (int i = 0; i < DiM; ++i)
    for (int j = 0; j < DiM; ++j)
      B[i][j] = A[i][j];
}

void CopyMat_withFactor(const MAT &A, MAT &B, double factor) {
  for (int i = 0; i < DiM; ++i)
    for (int j = 0; j < DiM; ++j)
      B[i][j] = factor * A[i][j];
}

void setValue(VEC &A, double value) {
  for (int i = 0; i < DiM; ++i)
    A[i] = value;
}

void setValue(MAT &A, double value) {
  for (int i = 0; i < DiM; ++i)
    for (int j = 0; j < DiM; ++j)
      A[i][j] = value;
}

void setValue(VECm1 &A, double value) {
#if DiM2
  A[0] = value;
#endif
#if DiM3
  A[0] = value;
  A[1] = value;
#endif
}

void setValue(MATm1 &A, double value) {
#if DiM2
  A[0][0] = value;
#endif
#if DiM3
  A[0][0] = value;
  A[0][1] = value;
  A[1][0] = value;
  A[1][1] = value;
#endif
}

void ProductMatVec(const MAT &matA, const VEC &vecB, VEC &matA_times_vecB) {
  double tmp;
  for (int i = 0; i < DiM; ++i) {
    tmp = 0.0;
    for (int j = 0; j < DiM; ++j)
      tmp += matA[i][j] * vecB[j];
    matA_times_vecB[i] = tmp;
  }
}

void ProductTransposeMatVec(const MAT &matA, const VEC &vecB,
                            VEC &matATranspose_times_vecB) {
  double tmp;
  for (int i = 0; i < DiM; ++i) {
    tmp = 0.0;
    for (int j = 0; j < DiM; ++j)
      tmp += matA[j][i] * vecB[j];
    matATranspose_times_vecB[i] = tmp;
  }
}

void ProductMatMat(const MAT &matA, const MAT &matB, MAT &matA_times_matB) {
  double tmp;
  for (int i = 0; i < DiM; ++i) {
    for (int j = 0; j < DiM; ++j) {
      tmp = 0.0;
      for (int k = 0; k < DiM; ++k)
        tmp += matA[i][k] * matB[k][j];
      matA_times_matB[i][j] = tmp;
    }
  }
}

void AddVec(const VEC &vecA, const VEC &vecB, VEC &vecA_plus_vecB) {
  for (int i = 0; i < DiM; ++i)
    vecA_plus_vecB[i] = vecA[i] + vecB[i];
}

void SubtractVec(const VEC &vecA, const VEC &vecB, VEC &vecA_minus_vecB) {
  for (int i = 0; i < DiM; ++i)
    vecA_minus_vecB[i] = vecA[i] - vecB[i];
}

void FactorVec(VEC &vec, double factor) {
  for (int i = 0; i < DiM; ++i)
    vec[i] *= factor;
}

void FactorMat(MAT &mat, double factor) {
  for (int i = 0; i < DiM; ++i)
    for (int j = 0; j < DiM; ++j)
      mat[i][j] *= factor;
}

void LinearCombination2Vecs(const VEC &vecA, const VEC &vecB, double factorA,
                            double factorB, VEC &fAvecAplusfBvecB) {
  for (int i = 0; i < DiM; ++i)
    fAvecAplusfBvecB[i] = factorA * vecA[i] + factorB * vecB[i];
}

double Norm2(const VEC &vec) {
  double nrm = vec[0] * vec[0];
#if DiM2a3_F
  nrm += vec[1] * vec[1];
#if DiM3
  nrm += vec[2] * vec[2];
#endif
#endif
  return sqrt(nrm);
}

double Normalize(const VEC &vec, VEC &normalized_vec) {
  double norm2 = Norm2(vec);
  double inv_norm = 1.0 / norm2;
  for (int i = 0; i < DiM; ++i)
    normalized_vec[i] = inv_norm * vec[i];
  return norm2;
}

#if DiM2a3_F
bool EigenValuesEigenVectors_PositiveSymMatrix(const MAT &mat,
                                               MAT &eigenVecsRowWise,
                                               VEC &eigenValues) {
#if DiM1
  eigenVecsRowWise[0][0] = 1.0;
  eigenValues[0] = mat[0][0];
  return true;
#endif
#if DiM2
  double a = mat[0][0], b = mat[1][1], c = mat[0][1];
  double tol = 1e-7 * (fabs(a) + fabs(b));
  double lambda0, lambda1;
  if (fabs(c) < tol) {
    if (a >= b) {
      lambda0 = a, lambda1 = b;
      eigenVecsRowWise[0][0] = 1.0;
      eigenVecsRowWise[0][1] = 0.0;
      eigenVecsRowWise[1][0] = 0.0;
      eigenVecsRowWise[1][1] = 1.0;
    } else {
      eigenVecsRowWise[0][0] = 0.0;
      eigenVecsRowWise[0][1] = 1.0;
      eigenVecsRowWise[1][0] = 1.0;
      eigenVecsRowWise[1][1] = 0.0;
    }
  } else {
    double amb = (a - b), apb = a + b;
    double sqrt_delta = sqrt(amb * amb + 4.0 * c * c);
    lambda0 = 0.5 * (apb + sqrt_delta),
    lambda1 = 0.5 * (apb - sqrt_delta); // decreasing order
    double theta = atan(c / (lambda0 - b));
    double costheta = cos(theta), sintheta = sin(theta);
    eigenVecsRowWise[0][0] = costheta;
    eigenVecsRowWise[0][1] = sintheta;
    eigenVecsRowWise[1][0] = -sintheta;
    eigenVecsRowWise[1][1] = costheta;
  }
  return true;
#endif
#if DiM3
  THROW("Eigenvalue functionis not implemented\n");
  return true;
#endif
}

double Norm2m1(const Vc_dm1 &vec) {
  double nrm = vec[0] * vec[0];
#if DiM3
  nrm += vec[1] * vec[1];
#endif
  return sqrt(nrm);
}

double Normalize_m1(const Vc_dm1 &vec, Vc_dm1 &normalized_vec) {
  double norm2 = Norm2m1(vec);
  double inv_norm = 1.0 / norm2;
  for (int i = 0; i < DiMm1; ++i)
    normalized_vec[i] = inv_norm * vec[i];
  return norm2;
}

void Extract_Shear(const VEC &vec_nt, Vc_dm1 &vec_t) {
  vec_t[0] = vec_nt[1];
#if DiM3
  vec_t[1] = vec_nt[2];
#endif
}
#endif

Vec_nt::Vec_nt() { nt_processed = false; }

void Vec_nt::SetZero() {
  setValue(vec, 0.0);
  nt_processed = false;
  vec_n = 0.0;
  vec_n_pos = 0.0;
#if DiM2a3_F
  setValue(vec_t, 0.0);
  vec_t_abs = 0.0;
  vec_effective = 0.0;
#endif
}

void Vec_nt::Compute_derivedValues_From_vec(const VEC &vecIn, double beta) {
  if (nt_processed)
    return;
  CopyVec(vecIn, vec);
  Compute_derivedValues_From_vec(beta);
  nt_processed = true;
}

void Vec_nt::Compute_derivedValues_From_vec(double beta) {
  vec_n = vec[0];
  vec_n_pos = 0.0;
  bool n_pos = (vec_n > 0.0);
  if (n_pos)
    vec_n_pos = vec_n;
#if DiM2a3_F
  Extract_Shear(vec, vec_t);
  vec_t_abs = Norm2m1(vec_t);
  vec_effective = vec_t_abs * beta;
  if (n_pos) {
    vec_effective = vec_effective * vec_effective;
    vec_effective += vec_n_pos * vec_n_pos;
    vec_effective = sqrt(vec_effective);
  }
#endif
}

void Vec_nt::Read(istream &in) {
  SetZero();
  string buf;
  READ_NSTRING(in, buf, buf);
  if (buf != "{") {
    if (buf == "}")
      return;
    else {
      THROW("istream should start with {");
    }
  }
  READ_NSTRING(in, buf, buf);
  while (buf != "}") {
    if (buf == "vec") {
      VEC vecIn;
      ReadV(vecIn, in);
      Compute_derivedValues_From_vec(vecIn, 1.0);
    } else if (buf == "nt_processed") {
      READ_NBOOL(in, buf, nt_processed);
    } else if (buf == "vec_n_pos") {
      READ_NDOUBLE(in, buf, vec_n_pos);
    }
#if DiM2a3_F
    else if (buf == "vec_t") {
      ReadVm1(vec_t, in);
    } else if (buf == "vec_t_abs") {
      READ_NDOUBLE(in, buf, vec_t_abs);
    } else if (buf == "vec_effective") {
      READ_NDOUBLE(in, buf, vec_effective);
    }
#endif
    else {
      cout << "buf:\t" << buf << '\n';
      THROW("invalid option\n");
    }
    READ_NSTRING(in, buf, buf);
  }
}
