#include "Configuration.h"

Configuration::Configuration() {
  wf_type = DG_2FUV;
  sOption = so_Riemann;
  dg_eps = dg_eps_p1;
  etaI = 1.0;
  etaB = 1.0;
  Initialize_Configuration();
}

void Configuration::Initialize_Configuration() {
  isDG = (wf_type != cfem_1D);
  weight_BLM_is_velocity = (wf_type == DG_2FUV);
  num_fields = 1;
  field_Start_U = 0;
  field_Start_V = 0;
  field_pos_weight_BLM = 0;
  if (wf_type == DG_2FUV) {
    num_fields = 2;
    field_pos_weight_BLM = 1;
  }
}

void Configuration::Compute_DG_Star_Weights_4_Inteior_Interface(
    const ElementProperties &left_ep, const ElementProperties &right_ep,
    bool insideDomainInterface, StarW2s &twoSideWeights) const {
  // calculates sigmaStar (not sigmaStar_n = sigmaStar . n) and wStar weight
  StarW1s shared_sigmaStar_wStarWeights;
  if (sOption == so_Riemann) {
    double Zl = left_ep.Z, Zr = right_ep.Z;
    double ZlpZr_inv = 1.0 / (Zl + Zr);
    double Zl_div_ZlpZr = Zl * Zl_div_ZlpZr, Zr_div_ZlpZr = Zr * Zl_div_ZlpZr;

    // sigma* weights
    shared_sigmaStar_wStarWeights.ss_f_sigmaL = Zl_div_ZlpZr;
    shared_sigmaStar_wStarWeights.ss_f_sigmaR = Zr_div_ZlpZr;
    shared_sigmaStar_wStarWeights.ss_f_wR = Zl * Zr * ZlpZr_inv;
    shared_sigmaStar_wStarWeights.ss_f_wL =
        -shared_sigmaStar_wStarWeights.ss_f_wR;

    // Velocity* weights (still not necessarily w weights)
    shared_sigmaStar_wStarWeights.ws_f_sigmaL = -ZlpZr_inv;
    shared_sigmaStar_wStarWeights.ws_f_sigmaR = ZlpZr_inv;
    shared_sigmaStar_wStarWeights.ws_f_wR = Zl_div_ZlpZr;
    shared_sigmaStar_wStarWeights.ws_f_wL = Zr_div_ZlpZr;
  } else {
    // taking care of the penality term
    if (etaI > 1e-15) {
      double eta = 0.0; // penality term - absolute value
      if (isW_velocity) {
        double Zl = left_ep.Z, Zr = right_ep.Z;
        double ZlpZr_inv = 1.0 / (Zl + Zr);
        eta = etaI * Zl * Zr * ZlpZr_inv;
      } else // w is u -> elliptic start
      {
        double ETilde = 0.5 * (left_ep.E + right_ep.E);
        double hTilde = 0.5 * (left_ep.hE + right_ep.hE);
        eta = etaI * ETilde / hTilde;
      }
      // sigma* = eta * (wR - wL)
      shared_sigmaStar_wStarWeights.ss_f_wR = -eta;
      shared_sigmaStar_wStarWeights.ss_f_wL = eta;
    }
    shared_sigmaStar_wStarWeights.ws_f_sigmaL = 0.0;
    shared_sigmaStar_wStarWeights.ws_f_sigmaR = 0.0;

    if (sOption == so_Central) {
      shared_sigmaStar_wStarWeights.ss_f_sigmaL = 0.5;
      shared_sigmaStar_wStarWeights.ss_f_sigmaR = 0.5;

      shared_sigmaStar_wStarWeights.ws_f_wR = 0.5;
      shared_sigmaStar_wStarWeights.ws_f_wL = 0.5;
    } else if (sOption == so_Alternating_sL) {
      shared_sigmaStar_wStarWeights.ss_f_sigmaL = 1.0;
      shared_sigmaStar_wStarWeights.ss_f_sigmaR = 0.0;

      shared_sigmaStar_wStarWeights.ws_f_wR = 0.0;
      shared_sigmaStar_wStarWeights.ws_f_wL = 1.0;
    } else if (sOption == so_Alternating_sR) {
      shared_sigmaStar_wStarWeights.ss_f_sigmaL = 0.0;
      shared_sigmaStar_wStarWeights.ss_f_sigmaR = 1.0;

      shared_sigmaStar_wStarWeights.ws_f_wR = 1.0;
      shared_sigmaStar_wStarWeights.ws_f_wL = 0.0;
    } else {
      cout << "sOption\t" << sOption << '\n';
      THROW("Invalid sOption\n");
    }
  }

  /// if w is u (Helmholtz), we need to use i omega to adjust the weights
  if (isHelmholtz) {
    Dcomplex iomega = omega * Icomp, iomega_inv = 1.0 / iomega;

    shared_sigmaStar_wStarWeights.ss_f_wR *= iomega;
    shared_sigmaStar_wStarWeights.ss_f_wL *= iomega;

    // Velocity* weights (still not necessarily w weights)
    shared_sigmaStar_wStarWeights.ws_f_sigmaL *= iomega_inv;
    shared_sigmaStar_wStarWeights.ws_f_sigmaR *= iomega_inv;
  }
  twoSideWeights.side_weights[SDL] = shared_sigmaStar_wStarWeights;
  twoSideWeights.side_weights[SDR] = shared_sigmaStar_wStarWeights;
  // need to multiply sigma*n parts of the right by -1, since for the right side
  // n = -1
  twoSideWeights.side_weights[SDR].ss_f_sigmaL =
      -twoSideWeights.side_weights[SDR].ss_f_sigmaL;
  twoSideWeights.side_weights[SDR].ss_f_sigmaR =
      -twoSideWeights.side_weights[SDR].ss_f_sigmaR;
  twoSideWeights.side_weights[SDR].ss_f_wL =
      -twoSideWeights.side_weights[SDR].ss_f_wL;
  twoSideWeights.side_weights[SDR].ss_f_wR =
      -twoSideWeights.side_weights[SDR].ss_f_wR;

  if ((!isBlochModeAnalysis) || insideDomainInterface) // ready to return
    return;

  // deal with Bloch mode analysis where gamma is involved
  twoSideWeights.side_weights[SDL].ss_f_sigmaR *= gamma;
  twoSideWeights.side_weights[SDL].ss_f_wR *= gamma;
  twoSideWeights.side_weights[SDL].ws_f_sigmaR *= gamma;
  twoSideWeights.side_weights[SDL].ws_f_wR *= gamma;

  twoSideWeights.side_weights[SDR].ss_f_sigmaL *= gamma_inv;
  twoSideWeights.side_weights[SDR].ss_f_wL *= gamma_inv;
  twoSideWeights.side_weights[SDR].ws_f_sigmaL *= gamma_inv;
  twoSideWeights.side_weights[SDR].ws_f_wL *= gamma_inv;
}
