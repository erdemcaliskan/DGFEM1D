#include "Configuration.h"
#include "commonMacros.h"
#include "float.h"
#include "globalFunctions.h"

#if USE_COMPLEX
#if !VCPP
#include <Eigen/Eigenvalues>
#endif
#endif

Configuration::Configuration()
{
	solutionMode = smt_Mats_OnlyNaturalMode;
    wf_type = DG_1F_vStar; // cfem DG_2FUV DG_1F_vStar;
    sOption = so_Riemann;
    dg_eps = dg_eps_p1;
	gamma = 1.0;
	gamma_inv = 1.0;
	omega = 1.0;

    etaI = 1.0;
    etaB = 1.0;
    num_elements = 5;
    xm = -0.5;
    xM = 0.5;
    E_uniform = 1.0;
    rho_uniform = 1.0;
    damping_uniform = 0.0;
    polyOrder = 1;
    leftBC = bct_Neumann;
    rightBC = bct_Neumann;
	leftBC_loadNumber = 0, rightBC_loadNumber = 0;
	leftBC_loadValues.push_back(0.0);
	rightBC_loadValues.push_back(1.0);

    outputFileNameWOExt = "output";
    lumpMass = false;

	isPeriodic = false;
	serialNumber = 0;
	mPropt = mpt_uniform;
	resolutionFactor = 0;
	sso_E = sso_mean_harmonic;
	sso_rho = sso_mean_arithmetic;
}

void Configuration::Main_SolveDomain(const string &configName, int serialNumberIn)
{
	serialNumber = serialNumberIn;
	serialNumber_str = "";
	if (serialNumber >= 0)
	{
		toString(serialNumber, serialNumber_str);
		serialNumber_str = "_" + serialNumber_str;
	}
    Read_ConfigFile(configName);
    Initialize_Configuration();
    Read_ElementGeometryProperties();
    Form_ElementMatrices();
    bool assembleMass = true; // if we solve a static probem we don't need this
    if (isDG)
        AssembleGlobalMatrices_DG(assembleMass);
    else
        AssembleGlobalMatrices_CFEM(assembleMass);
	Compute_StaticKaFSystem();
	Compute_BlochModeAnalysis();
    Process_Output_GlobalMatrices();
}

void Configuration::Initialize_Configuration()
{
	isDG = (wf_type != cfem);
	if (sz_omega > 0)
		solutionMode = smt_MatsF_Helmholtz;
	if (sz_k > 0)
		solutionMode = smt_Mats_OnlyBloch;
	if (solutionMode == smt_Mats_OnlyBloch)
	{
#if USE_COMPLEX
		if (isDG)
		{
			leftBC = bct_PeriodicOrBloch;
			rightBC = bct_PeriodicOrBloch;
		}
		else
		{
			leftBC = bct_Neumann;
			rightBC = bct_Neumann;
		}
#else
		THROW("solutionMode == smt_Mats_OnlyBloch is read as true from config file, but macro USE_COMPLEX is 0 in the code. In globalMacros.h turn it on\n");
#endif
	}
	b_PeriodicBC_But_NoBloch = (((leftBC == bct_PeriodicOrBloch) || (rightBC == bct_PeriodicOrBloch)) && (solutionMode != smt_Mats_OnlyBloch));

	SetBooleans_SolutionMode(solutionMode, sm_needForce, sm_isDynamic);

    weight_BLM_is_velocity = (wf_type == DG_2FUV);
    num_fields = 1;
    field_pos_weight_BLM = 0;
    if (wf_type == DG_2FUV)
    {
        num_fields = 2;
        field_pos_weight_BLM = 1;
    }
}

void Configuration::Compute_DG_Star_Weights_4_Inteior_Interface(const ElementProperties &left_ep,
                                                                const ElementProperties &right_ep,
                                                                bool insideDomainInterface,
                                                                StarW2s &twoSideWeights) const
{
    // calculates sigmaStar (not sigmaStar_n = sigmaStar . n) and wStar weight
    StarW1s shared_sigmaStar_wStarWeights;
    if (sOption == so_Riemann)
    {
        double Zl = left_ep.Z, Zr = right_ep.Z;
        double ZlpZr_inv = 1.0 / (Zl + Zr);
        double Zl_div_ZlpZr = Zl * ZlpZr_inv, Zr_div_ZlpZr = Zr * ZlpZr_inv;

        // sigma* weights
        shared_sigmaStar_wStarWeights.ss_f_sigmaL = Zr_div_ZlpZr;
        shared_sigmaStar_wStarWeights.ss_f_sigmaR = Zl_div_ZlpZr;
        shared_sigmaStar_wStarWeights.ss_f_wR = Zl * Zr * ZlpZr_inv;
        shared_sigmaStar_wStarWeights.ss_f_wL = -shared_sigmaStar_wStarWeights.ss_f_wR;

        // Velocity* weights (still not necessarily w weights) /// check this
        shared_sigmaStar_wStarWeights.ws_f_sigmaL = -ZlpZr_inv; // was + on 06/23/2023
        shared_sigmaStar_wStarWeights.ws_f_sigmaR = ZlpZr_inv;
        shared_sigmaStar_wStarWeights.ws_f_wR = Zr_div_ZlpZr;  // DAMPING
        shared_sigmaStar_wStarWeights.ws_f_wL = Zl_div_ZlpZr;  // DAMPING
    }
    else
    {
        // taking care of the penality term
        if (etaI > 1e-15)
        {
            double eta = 0.0; // penality term - absolute value
            if (IsWVelocity())
            { // w is velocity for 2F and 1Fv*
                double Zl = left_ep.Z, Zr = right_ep.Z;
                double ZlpZr_inv = 1.0 / (Zl + Zr);
                eta = etaI * Zl * Zr * ZlpZr_inv;
            }
            else // w is u -> elliptic start
            {
                double ETilde = 0.5 * (left_ep.E + right_ep.E);
                double hTilde = 0.5 * (left_ep.hE + right_ep.hE);
                eta = etaI * ETilde / hTilde;
            }
            // sigma* = eta * (wR - wL)
            shared_sigmaStar_wStarWeights.ss_f_wR =  eta; // was - 2023/06/23
            shared_sigmaStar_wStarWeights.ss_f_wL = -eta; // was + 2023/06/23
        }
        shared_sigmaStar_wStarWeights.ws_f_sigmaL = 0.0;
        shared_sigmaStar_wStarWeights.ws_f_sigmaR = 0.0;

        if (sOption == so_Central)
        {
            shared_sigmaStar_wStarWeights.ss_f_sigmaL = 0.5;
            shared_sigmaStar_wStarWeights.ss_f_sigmaR = 0.5;

            shared_sigmaStar_wStarWeights.ws_f_wR = 0.5;
            shared_sigmaStar_wStarWeights.ws_f_wL = 0.5;
        }
        else if (sOption == so_Alternating_sL)
        {
            shared_sigmaStar_wStarWeights.ss_f_sigmaL = 1.0;
            shared_sigmaStar_wStarWeights.ss_f_sigmaR = 0.0;

            shared_sigmaStar_wStarWeights.ws_f_wR = 0.0;
            shared_sigmaStar_wStarWeights.ws_f_wL = 1.0;
        }
        else if (sOption == so_Alternating_sR)
        {
            shared_sigmaStar_wStarWeights.ss_f_sigmaL = 0.0;
            shared_sigmaStar_wStarWeights.ss_f_sigmaR = 1.0;

            shared_sigmaStar_wStarWeights.ws_f_wR = 1.0;
            shared_sigmaStar_wStarWeights.ws_f_wL = 0.0;
        }
        else
        {
            cout << "sOption\t" << sOption << '\n';
            THROW("Invalid sOption\n");
        }
    }
	
    /// if w is u (Helmholtz), we need to use i omega to adjust the weights
#if USE_COMPLEX
    if (solutionMode == smt_MatsF_Helmholtz)
    {
        Dcomplex iomega = omega * Icomp, iomega_inv = 1.0 / iomega;

        shared_sigmaStar_wStarWeights.ss_f_wR *= iomega;
        shared_sigmaStar_wStarWeights.ss_f_wL *= iomega;

        // Velocity* weights (still not necessarily w weights)
        shared_sigmaStar_wStarWeights.ws_f_sigmaL *= iomega_inv;
        shared_sigmaStar_wStarWeights.ws_f_sigmaR *= iomega_inv;
    }
#endif
    twoSideWeights.side_weights[SDL] = shared_sigmaStar_wStarWeights;
    twoSideWeights.side_weights[SDR] = shared_sigmaStar_wStarWeights;
    // need to multiply sigma*n parts of the right by -1, since for the right side
    // n = -1
    twoSideWeights.side_weights[SDR].ss_f_sigmaL = -twoSideWeights.side_weights[SDR].ss_f_sigmaL;
    twoSideWeights.side_weights[SDR].ss_f_sigmaR = -twoSideWeights.side_weights[SDR].ss_f_sigmaR;
    twoSideWeights.side_weights[SDR].ss_f_wL = -twoSideWeights.side_weights[SDR].ss_f_wL;
    twoSideWeights.side_weights[SDR].ss_f_wR = -twoSideWeights.side_weights[SDR].ss_f_wR;

	
    if ((solutionMode != smt_Mats_OnlyBloch) || insideDomainInterface) // ready to return
        return;

    // deal with Bloch mode analysis where gamma is involved
#if USE_COMPLEX
    twoSideWeights.side_weights[SDL].ss_f_sigmaR *= gamma;
    twoSideWeights.side_weights[SDL].ss_f_wR *= gamma;
    twoSideWeights.side_weights[SDL].ws_f_sigmaR *= gamma;
    twoSideWeights.side_weights[SDL].ws_f_wR *= gamma;

    twoSideWeights.side_weights[SDR].ss_f_sigmaL *= gamma_inv;
    twoSideWeights.side_weights[SDR].ss_f_wL *= gamma_inv;
    twoSideWeights.side_weights[SDR].ws_f_sigmaL *= gamma_inv;
    twoSideWeights.side_weights[SDR].ws_f_wL *= gamma_inv;
#endif
}

bool Configuration::Compute_InteriorInterface_Matrices(const OneDimensionalElement &left_e,
                                                                        const OneDimensionalElement &right_e,
                                                                        bool insideDomainInterface,
                                                                        DCMATRIX &interfaceK,
                                                                        DCMATRIX &interfaceC) const
{
    StarW2s twoSideWeights;
    Compute_DG_Star_Weights_4_Inteior_Interface(left_e.elementProps, right_e.elementProps, insideDomainInterface,
                                                twoSideWeights);

    // sample numbers are provided for p = 2 -> ndof_parent_element = 3
    ////// We can move the calculation of start_ positions to the Initization step later

    unsigned int start_dof_l_usigma = 0;
    unsigned int start_dof_l_w; // w stands for u or v in the formulation
    unsigned int start_dof_r_usigma = 0;
    unsigned int start_dof_r_w;
    unsigned int size_interfaceMatrices = 2 * ndof_element;

    // matrix4wSlnField -> the stiffness terms that have velocity solution field ->
    //		for 1Fv* this will be interfaceC
    //		for other formulations 1Fu* and 2D this will be interfaceK
    interfaceK.resize(size_interfaceMatrices, size_interfaceMatrices);
    ZEROMAT(interfaceK, size_interfaceMatrices, size_interfaceMatrices);
    DCMATRIX *matrix4wSlnField;
    matrix4wSlnField = &interfaceK;
    bool hasC = false;

    double dg_eps_factor = (double)dg_eps;

    if (wf_type == DG_2FUV) // 2-field formulation
    {
        // ndof_element = 2 (u, v) * 3 (ndof_parent_element) = 6
        // size_interfaceMatrices = 2 (left and right) * ndof_element = 12
        // u and sigma of left element range of dof [0 -> 3]
        start_dof_l_usigma = 0;
        // v of left element range of dof [3 -> 6]
        start_dof_l_w = ndof_parent_element;
        // u and sigma of right element range of dof [6 -> 9]
        start_dof_r_usigma = 2 * ndof_parent_element;
        // v of right element range of dof [9 -> 12]
        start_dof_r_w = 3 * ndof_parent_element;
    }
    else
    {
        start_dof_l_usigma = 0;
        start_dof_l_w = 0;
        start_dof_r_usigma = ndof_parent_element;
        start_dof_r_w = start_dof_r_usigma;

        if (wf_type == DG_1F_vStar)
        {
            interfaceC.resize(size_interfaceMatrices, size_interfaceMatrices);
            ZEROMAT(interfaceC, size_interfaceMatrices, size_interfaceMatrices);
            matrix4wSlnField = &interfaceC;
            hasC = true;
            // time scale of the element is needed
            double beta_timeScale = 0.5 * (left_e.elementProps.time_e + right_e.elementProps.time_e);
            dg_eps_factor *= beta_timeScale;
        }
    }

    // N: shape: actual solution
    // B: dN/dx = B_xi * dx/dxi
    VECTOR N_leftElement(ndof_parent_element), N_rightElement(ndof_parent_element), B_leftElement(ndof_parent_element),
        B_rightElement(ndof_parent_element);
    double Jacobian_leftElement = left_e.GetJacobian();
    double Jacobian_rightElement = right_e.GetJacobian();
	double Jacobian_leftElementInv = 1.0 / Jacobian_leftElement;
	double Jacobian_rightElementInv = 1.0 / Jacobian_rightElement;

    for (unsigned int i = 0; i < ndof_parent_element; ++i)
    {
        N_leftElement[i] = parentElement.Ne_rightNode[i]; // * at the interface, shape of the left element is evaluated
                                                          // at its farthest right point -> Ne_rightNode
        N_rightElement[i] = parentElement.Ne_leftNode[i];

        B_leftElement[i] = Jacobian_leftElementInv * parentElement.Be_rightNode[i]; // similar to *
        B_rightElement[i] = Jacobian_rightElementInv * parentElement.Be_leftNode[i];
    }

    StarW1s *starFactors4LeftSide = &twoSideWeights.side_weights[SDL];
    StarW1s *starFactors4RightSide = &twoSideWeights.side_weights[SDR];
    // -wHat x sigma*_n - epsilon sigmaHat.n (w* - w)
    //! A: -wHat x sigma*_n | w = v for 2F and u for 1F formulations
    for (unsigned int i = 0; i < ndof_parent_element; ++i) // vHat of left element
    {
        /////////////////////////////// vHat from the left element
        // dVhatL_daL -> is the shape of the left element at its right point = N_leftElement
        double dVhatL_daL = N_leftElement[i];
        // dVhatR_daR -> is the shape of the right element at its left point = N_rigtElement
        double dVhatR_daR = N_rightElement[i];
        // a. Dependencies of Star on Left (sigmaL & velL)
        NUMBR dSigmaStarnL_dsigmaL = starFactors4LeftSide->ss_f_sigmaL;
        NUMBR dSigmaStarnL_dvelL = starFactors4LeftSide->ss_f_wL;
		NUMBR dSigmaStarnR_dsigmaL = starFactors4RightSide->ss_f_sigmaL;
		NUMBR dSigmaStarnR_dvelL = starFactors4RightSide->ss_f_wL;
		
		for (int unsigned j = 0; j < ndof_parent_element; ++j)
        {
            double d_sigmaL_daL = B_leftElement[j] * left_e.elementProps.E; // strain at that point times stress
            NUMBR dSigmaStar_daL = dSigmaStarnL_dsigmaL * d_sigmaL_daL;
            interfaceK(i + start_dof_l_w, j + start_dof_l_usigma) -= dVhatL_daL * dSigmaStar_daL;
			dSigmaStar_daL = dSigmaStarnR_dsigmaL * d_sigmaL_daL;
			interfaceK(i + start_dof_r_w, j + start_dof_l_usigma) -= dVhatR_daR * dSigmaStar_daL;

            double d_vell_daL = N_leftElement[j];
            dSigmaStar_daL = dSigmaStarnL_dvelL * d_vell_daL;
            (*matrix4wSlnField)(i + start_dof_l_w, j + start_dof_l_w) -= dVhatL_daL * dSigmaStar_daL;
			dSigmaStar_daL = dSigmaStarnR_dvelL * d_vell_daL;
			(*matrix4wSlnField)(i + start_dof_r_w, j + start_dof_l_w) -= dVhatR_daR * dSigmaStar_daL;
//            cout << interfaceK<< endl;
        }
        // b. Dependencies of Star on Right (sigmaR & velR)
        NUMBR dSigmaStarnL_dsigmaR = starFactors4LeftSide->ss_f_sigmaR;
        NUMBR dSigmaStarnL_dvelR = starFactors4LeftSide->ss_f_wR;
		NUMBR dSigmaStarnR_dsigmaR = starFactors4RightSide->ss_f_sigmaR;
		NUMBR dSigmaStarnR_dvelR = starFactors4RightSide->ss_f_wR;
		for (int unsigned j = 0; j < ndof_parent_element; ++j)
        {
            double d_sigmaR_daR = B_rightElement[j] * right_e.elementProps.E; // strain at that point times stress
            NUMBR dSigmaStar_daR = dSigmaStarnL_dsigmaR * d_sigmaR_daR;
            interfaceK(i + start_dof_l_w, j + start_dof_r_usigma) -= dVhatL_daL * dSigmaStar_daR;
			dSigmaStar_daR = dSigmaStarnR_dsigmaR * d_sigmaR_daR;
			interfaceK(i + start_dof_r_w, j + start_dof_r_usigma) -= dVhatR_daR * dSigmaStar_daR;

            double d_vell_daR = N_rightElement[j];
            dSigmaStar_daR = dSigmaStarnL_dvelR * d_vell_daR;
            (*matrix4wSlnField)(i + start_dof_l_w, j + start_dof_r_w) -= dVhatL_daL * dSigmaStar_daR;
			dSigmaStar_daR = dSigmaStarnR_dvelR * d_vell_daR;
			(*matrix4wSlnField)(i + start_dof_r_w, j + start_dof_r_w) -= dVhatR_daR * dSigmaStar_daR;
        }
    }
    DB(db << interfaceK << endl);
    //
    //! B: -epsilon (beta) sigmaHat.n (w* - w) [w = v for 2F, 1Fv, = u for 1Fu] -> beta = timeScale only for 1Fu*
    if (dg_eps == dg_eps_0)
        return hasC;

    for (unsigned int i = 0; i < ndof_parent_element; ++i) // vHat of left element
    {
        /////////////////////////////// sigmaHat from the left element
        double eps_x_d_sigmaL_dot_nL_daL =
            -dg_eps_factor * B_leftElement[i] * left_e.elementProps.E; // strain at that point times stress
        // nR = -1
        double eps_x_d_sigmaR_dot_nR_daR =
            dg_eps_factor * B_rightElement[i] * right_e.elementProps.E; // strain at that point times stress
        // a. dependence of w* - w = vel* - v(or u* - u) on the left element - [sigmaL & velL]
        NUMBR dwStarL_minus_w_dsigmaL = starFactors4LeftSide->ws_f_sigmaL;
        NUMBR dwStarL_minus_w_dwL = starFactors4LeftSide->ws_f_wL - 1.0; // (*) dw*/dwL - dwL/dwL = ws_f_wL - 1

		NUMBR dwStarR_minus_w_dsigmaL = starFactors4RightSide->ws_f_sigmaL;
		NUMBR dwStarR_minus_w_dwL = starFactors4RightSide->ws_f_wL; // (*) dw*/dwL - dwL/dwL = ws_f_wL - 1
		for (int unsigned j = 0; j < ndof_parent_element; ++j)
        {
            // sigma solutions
            double dSigmaL_daL = B_leftElement[j] * left_e.elementProps.E;
			// sigmaHatL.nL * (dw*/da)_left: dependency of w* on sigmaL
			NUMBR dwStar_minus_w_daL = dwStarL_minus_w_dsigmaL * dSigmaL_daL;
            interfaceK(i + start_dof_l_usigma, j + start_dof_l_usigma) +=
                eps_x_d_sigmaL_dot_nL_daL * dwStar_minus_w_daL;
            // sigmaHatR.nR * (dW*/da)_left: dependency of w* on sigmaL
			dwStar_minus_w_daL = dwStarR_minus_w_dsigmaL * dSigmaL_daL;
			interfaceK(i + start_dof_r_usigma, j + start_dof_l_usigma) +=
                eps_x_d_sigmaR_dot_nR_daR * dwStar_minus_w_daL;

            // w (u or v) solutions solutions
            double d_wL_daL = N_leftElement[j];
            dwStar_minus_w_daL = dwStarL_minus_w_dwL * d_wL_daL;
            // sigmaHatL.nL * (dw*/da)_left: dependency of w* on wL
            (*matrix4wSlnField)(i + start_dof_l_usigma, j + start_dof_l_w) +=
                eps_x_d_sigmaL_dot_nL_daL * dwStar_minus_w_daL;
            // sigmaHatR.nR * (dW*/da)_left: dependency of w* on wR
			dwStar_minus_w_daL = dwStarR_minus_w_dwL * d_wL_daL;
			(*matrix4wSlnField)(i + start_dof_r_usigma, j + start_dof_l_w) +=
                eps_x_d_sigmaR_dot_nR_daR * dwStar_minus_w_daL;
        }

        // b. dependence of w* = vel*(or u*) on the right element - [sigmaR & velR]
        NUMBR dwStarL_minus_w_dsigmaR = starFactors4LeftSide->ws_f_sigmaR;
        NUMBR dwStarL_minus_w_dwR = starFactors4LeftSide->ws_f_wR; // see (*) above

		NUMBR dwStarR_minus_w_dsigmaR = starFactors4RightSide->ws_f_sigmaR;
		NUMBR dwStarR_minus_w_dwR = starFactors4RightSide->ws_f_wR - 1.0; // see (*) above
		for (int unsigned j = 0; j < ndof_parent_element; ++j)
        {
            // sigma solutions
            double dSigmaR_daR = B_rightElement[j] * right_e.elementProps.E;
            NUMBR dwStar_minus_w_daR = dwStarL_minus_w_dsigmaR * dSigmaR_daR;
            // sigmaHatL.nL * (dw*/da)_right: dependency of w* on sigmaL
            interfaceK(i + start_dof_l_usigma, j + start_dof_r_usigma) +=
                eps_x_d_sigmaL_dot_nL_daL * dwStar_minus_w_daR;
            // sigmaHatR.nR * (dW*/da)_left: dependency of w* on sigmaL
			dwStar_minus_w_daR = dwStarR_minus_w_dsigmaR * dSigmaR_daR;
			interfaceK(i + start_dof_r_usigma, j + start_dof_r_usigma) +=
                eps_x_d_sigmaR_dot_nR_daR * dwStar_minus_w_daR;

            // w (u or v) solutions solutions
            double d_wR_daR = N_rightElement[j];
            dwStar_minus_w_daR = dwStarL_minus_w_dwR * d_wR_daR;
            // sigmaHatL.nL * (dw*/da)_left: dependency of w* on wL
            (*matrix4wSlnField)(i + start_dof_l_usigma, j + start_dof_r_w) +=
                eps_x_d_sigmaL_dot_nL_daL * dwStar_minus_w_daR;
            // sigmaHatR.nR * (dW*/da)_left: dependency of w* on wR
			dwStar_minus_w_daR = dwStarR_minus_w_dwR * d_wR_daR;
			(*matrix4wSlnField)(i + start_dof_r_usigma, j + start_dof_r_w) +=
                eps_x_d_sigmaR_dot_nR_daR * dwStar_minus_w_daR;
        }
    }
	return hasC;
}

bool Configuration::Compute_DirichletBoundaryInterface_Matrices(const OneDimensionalElement & e, bool leftDomainInterface, DCMATRIX & interfaceK, DCMATRIX & interfaceC) const
{
	// sample numbers are provided for p = 2 -> ndof_parent_element = 3
	unsigned int start_dof_usigma = 0;
	unsigned int start_dof_w; // w stands for u or v in the formulation
	unsigned int size_interfaceMatrices = ndof_element;

	// matrix4wSlnField -> the stiffness terms that have velocity solution field ->
	//		for 1Fv* this will be interfaceC
	//		for other formulations 1Fu* and 2D this will be interfaceK
	interfaceK.resize(size_interfaceMatrices, size_interfaceMatrices);
	ZEROMAT(interfaceK, size_interfaceMatrices, size_interfaceMatrices);
	DCMATRIX *matrix4wSlnField;
	matrix4wSlnField = &interfaceK;
	bool hasC = false;

	double dg_eps_factor = (double)dg_eps;
	if (wf_type == DG_2FUV) // 2-field formulation
		start_dof_w = ndof_parent_element;
	else
	{
		start_dof_w = 0;

		if (wf_type == DG_1F_vStar)
		{
			interfaceC.resize(size_interfaceMatrices, size_interfaceMatrices);
			ZEROMAT(interfaceC, size_interfaceMatrices, size_interfaceMatrices);
			matrix4wSlnField = &interfaceC;
			hasC = true;
			// time scale of the element is needed
			double beta_timeScale = e.elementProps.time_e;
			dg_eps_factor *= beta_timeScale;
		}
	}

	bool etaBZero = true;
	double etaFactor = 0.0;
	if (etaB > 1e-15)
	{
		if (IsWVelocity())
			etaFactor = etaB * 0.5 * e.elementProps.Z;
		else // w is u -> elliptic start
			etaFactor = etaB * e.elementProps.E/ e.elementProps.hE;
		etaBZero = false;
	}

	// N: shape: actual solution
	// B: dN/dx = B_xi * dx/dxi
	VECTOR N_lement(ndof_parent_element), B_Element(ndof_parent_element);
	double Jacobian_Element = e.GetJacobian();
	double n = 1.0;
	if (leftDomainInterface)
	{
		n = -1.0;
		for (unsigned int i = 0; i < ndof_parent_element; ++i)
		{
			N_lement[i] = parentElement.Ne_leftNode[i];
			B_Element[i] = parentElement.Be_leftNode[i] / Jacobian_Element;
		}
	}
	else
	{
		for (unsigned int i = 0; i < ndof_parent_element; ++i)
		{
			N_lement[i] = parentElement.Ne_rightNode[i];
			B_Element[i] = parentElement.Be_rightNode[i] / Jacobian_Element;
		}
	}

	// -wHat x sigma*_n - epsF * sigmaHat.n (wBar - w)

	//! A: -wHat x sigma*_n | w = v for 2F and u for 1F formulations
	//	sigma*n = sigma . n + etaFactor(wBar - w) ->
	//		A1: -wHat	x	sigma.n
	//		A2:	 wHat	x	etaFactor * w
	//! B: epsF * sigmaHat.n * w 
	for (unsigned int i = 0; i < ndof_parent_element; ++i)
	{
		double dVhat_da = N_lement[i];
		for (int unsigned j = 0; j < ndof_parent_element; ++j)
		{
			double d_sigma_da = B_Element[j] * e.elementProps.E; // strain at that point times stress
			NUMBR dSigmaStar_da = d_sigma_da * n; // sigma*n -> sigma . n  part
			NUMBR tmp = dVhat_da * dSigmaStar_da;
			// A1
			interfaceK(i + start_dof_w, j + start_dof_usigma) -= tmp;
			// B
			if (dg_eps != dg_eps_0)
				(*matrix4wSlnField)(j + start_dof_usigma, i + start_dof_w) += dg_eps_factor * tmp;

			// A2
			if (etaBZero == false)
				(*matrix4wSlnField)(i + start_dof_w, j + start_dof_w) += dVhat_da * etaFactor * N_lement[j];
		}
	}
	return hasC;
}

void Configuration::Compute_DirichletBoundaryInterface_ForceVector(const OneDimensionalElement &e,
	bool leftDomainInterface, double wBar, DCVECTOR& F) const
{
	// sample numbers are provided for p = 2 -> ndof_parent_element = 3
	unsigned int start_dof_usigma = 0;
	unsigned int start_dof_w; // w stands for u or v in the formulation
	unsigned int size_interfaceMatrices = ndof_element;

	ZEROVEC(F, size_interfaceMatrices);
	double dg_eps_factor = (double)dg_eps;
	if (wf_type == DG_2FUV) // 2-field formulation
		start_dof_w = ndof_parent_element;
	else
	{
		start_dof_w = 0;

		if (wf_type == DG_1F_vStar)
		{
			// time scale of the element is needed
			double beta_timeScale = e.elementProps.time_e;
			dg_eps_factor *= beta_timeScale;
		}
	}

	bool etaBZero = true;
	double etaFactor = 0.0;
	if (etaB > 1e-15)
	{
		if (IsWVelocity())
			etaFactor = etaB * 0.5 * e.elementProps.Z;
		else // w is u -> elliptic start
			etaFactor = etaB * e.elementProps.E / e.elementProps.hE;
		etaBZero = false;
	}

	// N: shape: actual solution
	// B: dN/dx = B_xi * dx/dxi
	VECTOR N_lement(ndof_parent_element), B_Element(ndof_parent_element);
	double Jacobian_Element = e.GetJacobian();
	double n = 1.0;
	if (leftDomainInterface)
	{
		n = -1.0;
		for (unsigned int i = 0; i < ndof_parent_element; ++i)
		{
			N_lement[i] = parentElement.Ne_leftNode[i];
			B_Element[i] = parentElement.Be_leftNode[i] / Jacobian_Element;
		}
	}
	else
	{
		for (unsigned int i = 0; i < ndof_parent_element; ++i)
		{
			N_lement[i] = parentElement.Ne_rightNode[i];
			B_Element[i] = parentElement.Be_rightNode[i] / Jacobian_Element;
		}
	}

	// -wHat x sigma*_n - epsF * sigmaHat.n (wBar - w)

	//! A: -wHat x sigma*_n | w = v for 2F and u for 1F formulations
	//	sigma*n = sigma . n + etaFactor(wBar - w)			-> No Force
	//		A1: -wHat	x	sigma.n							
	//		A2:	-wHat	x	etaFactor * (wBar - w)			-> wHat	x	etaFactor * wBar 
	//! B: -epsF * sigmaHat.n * (wBar - w)					-> epsF * sigmaHat.n *  wBar 
	// B
	if (dg_eps != dg_eps_0)
	{
		for (unsigned int i = 0; i < ndof_parent_element; ++i) // vHat of left element
		{
			double d_sigma_da = B_Element[i] * e.elementProps.E; // strain at that point times stress
			double dSigmaStar_da = d_sigma_da * n; // sigma*n -> sigma . n  part
			F(i + start_dof_usigma) += dg_eps_factor * dSigmaStar_da * wBar;
		}
	}
	// A2
	if (etaBZero == false)
	{
		for (unsigned int i = 0; i < ndof_parent_element; ++i) // vHat of left element
		{
			double dVhat_da = N_lement[i];
			F(i + start_dof_w) += dVhat_da * etaFactor * wBar;
		}
	}
}

void Configuration::Compute_NeumannBoundaryInterface_ForceVector(const OneDimensionalElement &e,
	bool leftDomainInterface, double sigmanBarn, DCVECTOR& F) const
{
	// sample numbers are provided for p = 2 -> ndof_parent_element = 3
	unsigned int start_dof_usigma = 0;
	unsigned int start_dof_w = 0; // w stands for u or v in the formulation
	if (wf_type == DG_2FUV) // 2-field formulation
		start_dof_w = ndof_parent_element;

	unsigned int size_interfaceMatrices = ndof_element;
	ZEROVEC(F, size_interfaceMatrices);

	// N: shape: actual solution, B not needed
	VECTOR N_lement(ndof_parent_element);
	double Jacobian_Element = e.GetJacobian();
	//double n = 1.0;
	if (leftDomainInterface)
	{
		//n = -1.0;
		for (unsigned int i = 0; i < ndof_parent_element; ++i)
			N_lement[i] = parentElement.Ne_leftNode[i];
	}
	else
	{
		for (unsigned int i = 0; i < ndof_parent_element; ++i)
			N_lement[i] = parentElement.Ne_rightNode[i];
	}

	// wHat x sigmaBarn
	for (unsigned int i = 0; i < ndof_parent_element; ++i) // vHat of left element
		F(i + start_dof_w) = N_lement[i] * sigmanBarn;
}

void Configuration::Compute_left_right_prescribedVals(double time, double & valLeft, double & valRight)
{
	valLeft = ComputeLoad(time, leftBC_loadNumber, leftBC_loadValues);
	valRight = ComputeLoad(time, rightBC_loadNumber, rightBC_loadValues);
}

bool Configuration::HasNonZeroDampingMatrix() const
{
    if (wf_type == DG_2FUV)
        return false;
    if (wf_type == DG_1F_vStar)
        return true;

    static double tol_4_damping = 1000 * DBL_MIN;
    if ((wf_type == cfem) || (wf_type == DG_1F_uStar))
    {
        // depends on the material damping
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            if (elements[ei].elementProps.damping > tol_4_damping)
                return true;
        }
        return false;
    }
    cout << "wf_type\t" << wf_type << '\n';
    THROW("wf_type case not implemented\n");
}

void Configuration::Read_ElementGeometryProperties()
{
    // for simple case only needed
    if (mPropt == mpt_uniform)
    {
		if (num_elements <= 0)
		{
			THROW("(num_elements <= 0)\n");
		}
		domain_length = (xM - xm);
        double hE = domain_length / num_elements;
        elements.resize(num_elements);
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            elements[ei].elementProps.hE = hE;
            elements[ei].elementProps.E  = E_uniform;
            elements[ei].elementProps.rho = rho_uniform;
            elements[ei].elementProps.damping = damping_uniform;
        }
    }
    else if (mPropt == mpt_random_field)
    {
		bool hasRandom_E = !randomField_E.IsEmpty(), hasRandom_rho = !randomField_rho.IsEmpty();
		vector<double> xs;
		if (hasRandom_E)
		{
			num_elements = randomField_E.getNumSegments();
			randomField_E.get_xs(xs);
		}
		else if (hasRandom_rho)
		{
			num_elements = randomField_rho.getNumSegments();
			randomField_rho.get_xs(xs);
		}
		else
		{
			THROW("Neither E or rho is read from a field and cannot set the domain properties for (mPropt == mpt_random_field)\n");
		}
		elements.resize(num_elements);
		bool has_r_E = ((hasRandom_E) && (sso_E != sso_none));
		bool has_r_rho = ((hasRandom_rho) && (sso_rho != sso_none));
		for (unsigned int ei = 0; ei < num_elements; ++ei)
		{
			elements[ei].elementProps.hE = xs[ei + 1] - xs[ei];
			if (has_r_E)
				elements[ei].elementProps.E = randomField_E.getSegmentValueByIndex(ei);
			else
				elements[ei].elementProps.E = E_uniform;

			if (has_r_rho)
				elements[ei].elementProps.rho = randomField_rho.getSegmentValueByIndex(ei);
			else
				elements[ei].elementProps.rho = rho_uniform;

			elements[ei].elementProps.damping = damping_uniform;
		}
    }
	else if (mPropt == mpt_layered)
	{
		num_elements = layered_properties.sz_allsequences;
		elements.resize(num_elements);
		xm = 0.0, xM = 0.0;
		double h;
		for (unsigned int ei = 0; ei < num_elements; ++ei)
		{
			oneBulk_Elastic_Prop* bepPtr = &layered_properties.finalBulkSegments[ei];
			h = bepPtr->getLength();
			xM += h;
			elements[ei].elementProps.hE = h;
			elements[ei].elementProps.E = bepPtr->E;
			elements[ei].elementProps.rho = bepPtr->rho;
			elements[ei].elementProps.damping = bepPtr->damping;
		}
		double halfLength = 0.5 * xM;
		xm -= halfLength;
		xM -= halfLength;
	}
	else
	{
		cout << "mPropt\t" << mPropt << '\n';
		THROW("Invalid mPropt\n");
	}
	domain_length = xM - xm;
}

void Configuration::Form_ElementMatrices()
{
    // initialize the parent
    parentElement.Initialize(polyOrder, lumpMass);
    ndof_parent_element = parentElement.ndof;
    ndof_element = ndof_parent_element * num_fields;

    edof_2_globalDofMap.resize(num_elements);
    element_start_dof.resize(num_elements);
    element_end_dof.resize(num_elements);

    edof_2_globalDofMap.clear();

    if (wf_type != cfem)
    {
        nfdof_domain = ndof_element * num_elements;
        ndof_domain = nfdof_domain;
        npdof_domain = 0;
        vector<int> dof_map_element(ndof_element);
        unsigned cntr = 0;
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            element_start_dof[ei] = cntr;
            element_end_dof[ei] = cntr + ndof_element;
            for (unsigned int j = 0; j < ndof_element; ++j)
                dof_map_element[j] = cntr++;
            edof_2_globalDofMap.push_back(dof_map_element);
        }
    }
    else
    {
        /// A. Number of dofs and domain:nodal dof map
        // there is one dof shared between elements
        number_of_nodes = (ndof_element - 1) * num_elements + 1;
        //
        nodedof_map.resize(number_of_nodes);
        // for non-periodic case, there is one extra dof on far right
        if (!b_PeriodicBC_But_NoBloch)
            ndof_domain = number_of_nodes;
        else
            ndof_domain = number_of_nodes - 1;

        int cntr_f = 0, cntr_p = 0;
        if (leftBC == bct_Dirichlet)
            nodedof_map[0] = -(++cntr_p); // it will make it 1 and set -1 to dof of node 0
        else
            nodedof_map[0] = cntr_f++;

        if (rightBC == bct_Dirichlet)
            nodedof_map[ndof_domain - 1] = -(++cntr_p); // it will make it 1 and set -1 to dof of node 0

        // how many prescribed dofs we have
        npdof_domain = cntr_p;
        nfdof_domain = ndof_domain - npdof_domain;

        unsigned int en = ndof_domain - 2 + (int)b_PeriodicBC_But_NoBloch;
        // all of these are free
        for (unsigned int nodei = 1; nodei <= en; ++nodei)
            nodedof_map[nodei] = cntr_f++;
        // finally the far right dof/node
        if (rightBC != bct_Dirichlet)
        {
            unsigned int lastNode_index = number_of_nodes - 1;
            if (b_PeriodicBC_But_NoBloch)
                nodedof_map[lastNode_index] = 0;
            else if ((rightBC == bct_Neumann) ||
                     (rightBC == bct_PeriodicOrBloch)) // for the latter it will be Bloch as periodic is handled above.
                                                       // Bloch is handled as Neuman
                nodedof_map[lastNode_index] = cntr_f++;
            else
            {
                cout << "rightBC\t" << rightBC << '\n';
                THROW("Right BC is not taken care of\n");
            }
        }

        /// B. Element dof maps
        unsigned int st;
        edof_2_globalDofMap.resize(num_elements);
        int ndof_elementm1 = (ndof_element - 1); // number of "unshared" dofs per element
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            st = ei * ndof_elementm1; 
            vector<int> *e_dofMap = &edof_2_globalDofMap[ei];
            e_dofMap->resize(ndof_element);
            (*e_dofMap)[0] = nodedof_map[st]; // left node
            (*e_dofMap)[1] = nodedof_map[st + ndof_element - 1]; // left node
 
            for (unsigned int j = 2; j < ndof_element; ++j)
                (*e_dofMap)[j] = nodedof_map[st + j - 1];
        }
    }

    b_hasNonZeroDampingMatrix = HasNonZeroDampingMatrix();

    for (unsigned int ei = 0; ei < num_elements; ++ei)
    {
        OneDimensionalElement *ePtr = &elements[ei];
        ePtr->elementProps.Initialize_ElementProperties();

        // these building blocks are mpe, kpe with 1F size but with real element geometry (size) and no material
        // property
        MATRIX meBuildingBlock, keBuildingBlock;
        double Jacobian = 0.5 * ePtr->elementProps.hE;
        double JacobianInv = 1.0 / Jacobian; 
        keBuildingBlock.Multiply(parentElement.kpe, JacobianInv);
        meBuildingBlock.Multiply(parentElement.mpe, Jacobian);

        ePtr->ke.resize(ndof_element, ndof_element);
        ePtr->ke = 0.0;
        ePtr->me.resize(ndof_element, ndof_element);
        ePtr->me = 0.0;
        if (b_hasNonZeroDampingMatrix)
        {
            ePtr->ce.resize(ndof_element, ndof_element);
            ePtr->ce = 0.0;
        }

        if (wf_type == DG_2FUV)
        {
            double alpha = 1.0 / ePtr->elementProps.time_e;
            alpha *= alpha;
            for (unsigned int i = 0; i < parentElement.ndof; ++i)
            {
                for (unsigned int j = 0; j < parentElement.ndof; ++j)
                {
                    /// stiffness term
                    // -alpha uHat * v
                    ePtr->ke[i][j + parentElement.ndof] = -alpha * meBuildingBlock[i][j];
                    // vHat * damping * v
                    ePtr->ke[i + parentElement.ndof][j + parentElement.ndof] =
                        ePtr->elementProps.damping * meBuildingBlock[i][j];
                    // grad vHat * sigma
                    ePtr->ke[i + parentElement.ndof][j] = ePtr->elementProps.E * keBuildingBlock[i][j];

                    /// mass term
                    // alpha uHat * uDot
                    ePtr->me[i][j] = alpha * meBuildingBlock[i][j];
                    // vHat * vDot
                    ePtr->me[i + parentElement.ndof][j + parentElement.ndof] =
                        ePtr->elementProps.rho * meBuildingBlock[i][j];
                }
            }
        }
        else // if ((wf_type == DG_1F_vStar) || (wf_type == DG_1F_uStar) || (wf_type == cfem))
        {
            for (unsigned int i = 0; i < parentElement.ndof; ++i)
            {
                for (unsigned int j = 0; j < parentElement.ndof; ++j)
                {
                    /// stiffness term
                    // grad uHat * sigma
                    ePtr->ke[i][j] = ePtr->elementProps.E * keBuildingBlock[i][j];
                    /// mass term
                    // uHat * vDot = uHat * uDDot
                    ePtr->me[i][j] = ePtr->elementProps.rho * meBuildingBlock[i][j];

                    /// damping term
                    // uHat * damping * v = uHat * damping * uDot
                    if (b_hasNonZeroDampingMatrix)
                        ePtr->ce[i][j] = ePtr->elementProps.damping * meBuildingBlock[i][j];
                }
            }
        }
    }
}

void Configuration::AssembleGlobalMatrices_DG(bool assembleMassIn)
{
    globalK.resize(nfdof_domain, nfdof_domain);
    
	ZEROMAT(globalK,nfdof_domain,nfdof_domain);

	if (assembleMassIn)
    {
        globalM.resize(nfdof_domain, nfdof_domain);
		ZEROMAT(globalM,nfdof_domain,nfdof_domain);
    }
    if (b_hasNonZeroDampingMatrix)
    {
        globalC.resize(nfdof_domain, nfdof_domain);
		ZEROMAT(globalC,nfdof_domain,nfdof_domain);
    }

    /// Step A:  Interior of the elments
    unsigned int st;
    for (unsigned int ei = 0; ei < num_elements; ++ei)
    {
        st = element_start_dof[ei];
        OneDimensionalElement *ePtr = &elements[ei];
        for (unsigned int i = 0; i < ndof_element; ++i)
            for (unsigned int j = 0; j < ndof_element; ++j)
                globalK(i + st, j + st) += ePtr->ke[i][j];
    }
    if (assembleMassIn)
    {
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            st = element_start_dof[ei];
            OneDimensionalElement *ePtr = &elements[ei];
            for (unsigned int i = 0; i < ndof_element; ++i)
                for (unsigned int j = 0; j < ndof_element; ++j)
                    globalM(i + st, j + st) += ePtr->me[i][j];
        }
    }
    if (b_hasNonZeroDampingMatrix)
    {
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            st = element_start_dof[ei];
            OneDimensionalElement *ePtr = &elements[ei];
            for (unsigned int i = 0; i < ndof_element; ++i)
                for (unsigned int j = 0; j < ndof_element; ++j)
                {
                    globalC(i + st, j + st) += ePtr->ce[i][j];
                }
        }
    }
//    }

    /// Step B:  Inter-element facets

    if (leftBC == bct_Characteristics)
        THROW("left BC has stiffness contributions not taken care of\n");
    if (rightBC == bct_Characteristics)
        THROW("right BC has stiffness contributions not taken care of\n");

	unsigned int en = num_elements - 1;

	// Dirichlet & Neumman BC Initialization and  K / C
	e_DirichletBCs.clear();
	leftDomainInterface_DirichletBCs.clear();
	eIndex_DirichletBCs.clear();

	e_NeumannBCs.clear();
	leftDomainInterface_NeumannBCs.clear();
	eIndex_NeumannBCs.clear();

	if (leftBC == bct_Dirichlet)
	{
		e_DirichletBCs.push_back(&elements[0]);
		leftDomainInterface_DirichletBCs.push_back(true);
		eIndex_DirichletBCs.push_back(0);
	}
	else if (leftBC == bct_Neumann)
	{
		e_NeumannBCs.push_back(&elements[0]);
		leftDomainInterface_NeumannBCs.push_back(true);
		eIndex_NeumannBCs.push_back(0);
	}
	if (rightBC == bct_Dirichlet)
	{
		e_DirichletBCs.push_back(&elements[en]);
		leftDomainInterface_DirichletBCs.push_back(false);
		eIndex_DirichletBCs.push_back(en);
	}
	else if (rightBC == bct_Neumann)
	{
		e_NeumannBCs.push_back(&elements[en]);
		leftDomainInterface_NeumannBCs.push_back(false);
		eIndex_NeumannBCs.push_back(en);
	}
	num_DirichletBCs = e_DirichletBCs.size();
	num_NeumannBCs = e_NeumannBCs.size();
	for (unsigned int dbci = 0; dbci < num_DirichletBCs; ++dbci)
	{
		DCMATRIX interfaceK, interfaceC;
		ZEROMAT(interfaceK, ndof_element, ndof_element);
		ZEROMAT(interfaceC, ndof_element, ndof_element);
		bool hasC = Compute_DirichletBoundaryInterface_Matrices(*(e_DirichletBCs[dbci]), leftDomainInterface_DirichletBCs[dbci], interfaceK, interfaceC);

		unsigned int pos = eIndex_DirichletBCs[dbci];
		unsigned int st = element_start_dof[pos];

		for (unsigned int i = 0; i < ndof_element; ++i)
		{
			for (unsigned int j = 0; j < ndof_element; ++j)
				globalK(i + st, j + st) += interfaceK(i, j);
		}
		if (hasC)
		{
			for (unsigned int i = 0; i < ndof_element; ++i)
			{
				for (unsigned int j = 0; j < ndof_element; ++j)
					globalC(i + st, j + st) += interfaceC(i, j);
			}
		}
	}

    // n: number of elements -> n + 1 interfaces
    // I0 eInterior0 I1 eInterior1 ... I(n-1) eInterior(n-1) In
    // interior interfaces are from 1 to n -1
    st = 1;
	unsigned enInterface = en;
    if (b_PeriodicBC_But_NoBloch)
		enInterface = num_elements;

	unsigned int interface_dof = 2 * ndof_element;

    for (unsigned int interfacei = st; interfacei <= enInterface; ++interfacei)
    {
        const OneDimensionalElement *left_ePtr, *right_ePtr;
        int left_e_dof_st, right_e_dof_st;
        bool insideDomainInterface = true;
		left_e_dof_st = element_start_dof[interfacei - 1];
		left_ePtr = &elements[interfacei - 1];
		if (interfacei < num_elements) // interior interfaces
        {
			right_e_dof_st = element_start_dof[interfacei];
			right_ePtr = &elements[interfacei];
		}
        else
        {
			right_e_dof_st = element_start_dof[0];
			right_ePtr = &elements[0];
			insideDomainInterface = false;
        }
        DCMATRIX interfaceK, interfaceC;
		ZEROMAT(interfaceK, interface_dof, interface_dof);
		ZEROMAT(interfaceC, interface_dof, interface_dof);

        bool hasC = Compute_InteriorInterface_Matrices(*left_ePtr, *right_ePtr, insideDomainInterface,
                                                                        interfaceK, interfaceC);

        /// assemble K
        if (insideDomainInterface) // interior interfaces
        {
            for (unsigned int i = 0; i < interface_dof; ++i)
            {
                for (unsigned int j = 0; j < interface_dof; ++j)
                {
                    globalK(i + left_e_dof_st, j + left_e_dof_st) += interfaceK(i, j);
//                    cout << "K(" << i + left_e_dof_st << "," << j + left_e_dof_st << ") = " << interfaceK(i, j) << endl;
                }
            }
            if (hasC)
            {
                for (unsigned int i = 0; i < interface_dof; ++i)
                {
                    for (unsigned int j = 0; j < interface_dof; ++j)
                    {
                        globalC(i + left_e_dof_st, j + left_e_dof_st) += interfaceC(i, j);
//                        cout << "interfaceC = " << interfaceC << endl;
//                        cout << "C = " << globalC << endl;
                    }
                }
            }
        }
        else // This part is only done for period (but not Bloch) BC. Bloch BC is done later for all ks
        {
            for (unsigned int i = 0; i < ndof_element; ++i)
            {
                for (unsigned int j = 0; j < ndof_element; ++j)
                {
                    // RR
                    globalK(i + right_e_dof_st, j + right_e_dof_st) += interfaceK(i + ndof_element, j + ndof_element);
                    // RL
                    globalK(i + right_e_dof_st, j + left_e_dof_st) += interfaceK(i + ndof_element, j);
                    // LR
                    globalK(i + left_e_dof_st, j + right_e_dof_st) += interfaceK(i, j + ndof_element);
                    // LL
                    globalK(i + left_e_dof_st, j + left_e_dof_st) += interfaceK(i, j);
				}
            }
            if (hasC)
            {
                for (unsigned int i = 0; i < ndof_element; ++i)
                {
                    for (unsigned int j = 0; j < ndof_element; ++j)
                    {
                        // RR
						globalC(i + right_e_dof_st, j + right_e_dof_st) += interfaceC(i + ndof_element, j + ndof_element);
						// RL
						globalC(i + right_e_dof_st, j + left_e_dof_st) += interfaceC(i + ndof_element, j);
						// LR
						globalC(i + left_e_dof_st, j + right_e_dof_st) += interfaceC(i, j + ndof_element);
						// LL
						globalC(i + left_e_dof_st, j + left_e_dof_st) += interfaceC(i, j);
					}
                }
            }
        }
    }
}

void Configuration::AssembleGlobalMatrices_CFEM(bool assembleMassIn)
{
    globalK.resize(nfdof_domain, nfdof_domain);
	ZEROMAT(globalK,nfdof_domain,nfdof_domain);

    if (assembleMassIn)
    {
        globalM.resize(nfdof_domain, nfdof_domain);
		ZEROMAT(globalM,nfdof_domain,nfdof_domain);
    }
    if (b_hasNonZeroDampingMatrix)
    {
        globalC.resize(nfdof_domain, nfdof_domain);
		ZEROMAT(globalC,nfdof_domain,nfdof_domain);
    }

    int I, J; // I, J are global matrix rows and columns
    // i, j are local matrix rows and columns

    /// Step A:  Interior of the elments
    for (unsigned int ei = 0; ei < num_elements; ++ei)
    {
        OneDimensionalElement *ePtr = &elements[ei];
        vector<int> *e_dofMap = &edof_2_globalDofMap[ei];
        for (unsigned int i = 0; i < ndof_element; ++i)
        {
            I = (*e_dofMap)[i];
            for (unsigned int j = 0; j < ndof_element; ++j)
            {
                J = (*e_dofMap)[j];
                if ( double(I) >= 0.0 && double(J) >= 0 )
                {
                    globalK(I, J) += ePtr->ke[i][j];
                }
 //               globalK(I, J) += ePtr->ke[i][j];
            }
        }
    }
    if (assembleMassIn)
    {
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            OneDimensionalElement *ePtr = &elements[ei];
            vector<int> *e_dofMap = &edof_2_globalDofMap[ei];
            for (unsigned int i = 0; i < ndof_element; ++i)
            {
                I = (*e_dofMap)[i];
                for (unsigned int j = 0; j < ndof_element; ++j)
                {
                    J = (*e_dofMap)[j];
                    if ( I >= 0 && J >= 0 )
                        globalM(I, J) += ePtr->me[i][j];
                }
            }
        }
    }

    if (b_hasNonZeroDampingMatrix)
    {
        for (unsigned int ei = 0; ei < num_elements; ++ei)
        {
            OneDimensionalElement *ePtr = &elements[ei];
            vector<int> *e_dofMap = &edof_2_globalDofMap[ei];
            for (unsigned int i = 0; i < ndof_element; ++i)
            {
                I = (*e_dofMap)[i];
                for (unsigned int j = 0; j < ndof_element; ++j)
                {
                    J = (*e_dofMap)[j];
                    globalC(I, J) += ePtr->ce[i][j];
                }
            }
        }
    }
}

void Configuration::Read_ConfigFile(const string& configName)
{
    if (configName == "")
        return;
    string buf;
    fstream in(configName.c_str(), ios::in);
    if (!in.is_open())
    {
        cout << "configName\t" << configName << '\n';
        THROW("cannot open config file\n");
    }
    READ_NSTRING(in, buf, buf);
    if (buf != "{")
    {
        if (buf == "}")
            return;
        else
        {
            THROW("istream should start with {");
        }
    }
	VecOfVals ks_real, ks_imag, omegas_real, omegas_imag;

    READ_NSTRING(in, buf, buf);
    while (buf != "}")
    {
        if (buf == "polyOrder")
        {
            READ_NINTEGER(in, buf, polyOrder);
        }
        else if (buf == "lumpMass")
        {
            READ_NBOOL(in, buf, lumpMass);
        }
		else if (buf == "solutionMode")
		{
			in >> solutionMode;
		}
		else if (buf == "wf_type")
		{
			in >> wf_type;
		}
		else if (buf == "ks_real")
		{
			in >> ks_real;
		}
		else if (buf == "ks_imag")
		{
			in >> ks_imag;
		}
		else if (buf == "omegas_real")
		{
			in >> omegas_real;
		}
		else if (buf == "omegas_imag")
		{
			in >> omegas_imag;
		}
		else if (buf == "sOption")
		{
			in >> sOption;
		}
		else if (buf == "dg_eps")
		{
			int tmpi;
			READ_NINTEGER(in, buf, tmpi);
			dg_eps = DG_epsilonT(tmpi);
		}
		else if (buf == "etaI")
		{
			READ_NDOUBLE(in, buf, etaI);
		}
		else if (buf == "etaB")
		{
			READ_NDOUBLE(in, buf, etaB);
		}
		else if (buf == "leftBC")
		{
			in >> leftBC;
		}
		else if (buf == "rightBC")
		{
			in >> rightBC;
		}
		else if (buf == "leftBC_loadNumber")
		{
			READ_NINTEGER(in, buf, leftBC_loadNumber);
		}
		else if (buf == "rightBC_loadNumber")
		{
			READ_NINTEGER(in, buf, rightBC_loadNumber);
		}
		else if (buf == "leftBC_loadValues")
		{
			ReadVectorDouble(in, leftBC_loadValues);
		}
		else if (buf == "rightBC_loadValues")
		{
			ReadVectorDouble(in, rightBC_loadValues);
		}
		else if (buf == "outputFileNameWOExt")
        {
            READ_NSTRING(in, buf, outputFileNameWOExt);
        }
		else if (buf == "mPropt")
		{
			in >> mPropt;
		}
		else if (buf == "num_elements")
		{
			READ_NINTEGER(in, buf, num_elements);
		}
		else if (buf == "xm")
		{
			READ_NDOUBLE(in, buf, xm);
		}
		else if (buf == "xM")
		{
			READ_NDOUBLE(in, buf, xM);
		}
		else if (buf == "E_uniform")
		{
			READ_NDOUBLE(in, buf, E_uniform);
		}
		else if (buf == "rho_uniform")
		{
			READ_NDOUBLE(in, buf, rho_uniform);
		}
		else if (buf == "damping_uniform")
		{
			READ_NDOUBLE(in, buf, damping_uniform);
		}
		else if (buf == "resolutionFactor")
		{
			READ_NINTEGER(in, buf, resolutionFactor);
		}
		else if (buf == "sso_E")
		{
			in >> sso_E;
		}
		else if (buf == "sso_rho")
		{
			in >> sso_rho;
		}
		else if (buf == "sso_rho")
		{
			in >> sso_rho;
		}
		else if (buf == "randomField_E")
		{
			string baseName_WOExt;
			READ_NSTRING(in, buf, baseName_WOExt);
			string dataFileName = baseName_WOExt + serialNumber_str + ".txt";
			fstream inData;
			bool readData = ((mPropt == mpt_random_field) && (sso_E != sso_none));
			if (readData)
			{
				inData.open(dataFileName.c_str(), ios::in);
				if (!inData.is_open())
				{
					cout << "dataFileName\t" << dataFileName << '\n';
					THROW("Cannot open file\n");
				}
			}
			randomField_E.Read_Initialize_OneIHField(inData, &in, &isPeriodic, &xM, &xm, resolutionFactor, sso_E, readData);
		}
		else if (buf == "randomField_rho")
		{
			string baseName_WOExt;
			READ_NSTRING(in, buf, baseName_WOExt);
			string dataFileName = baseName_WOExt + serialNumber_str + ".txt";
			fstream inData;
			bool readData = ((mPropt == mpt_random_field) && (sso_rho != sso_none));
			if (readData)
			{
				inData.open(dataFileName.c_str(), ios::in);
				if (!inData.is_open())
				{
					cout << "dataFileName\t" << dataFileName << '\n';
					THROW("Cannot open file\n");
				}
			}
			randomField_rho.Read_Initialize_OneIHField(inData, &in, &isPeriodic, &xM, &xm, resolutionFactor, sso_rho, readData);
		}
		else if (buf == "layered_properties")
		{
			bool readData = (mPropt == mpt_layered);
			layered_properties.Read_Bulk_Elastic_Prop(in, serialNumber, readData);
		}
        else
        {
            cout << "buf:\t" << buf << '\n';
            THROW("invalid option\n");
        }
        READ_NSTRING(in, buf, buf);
    }
	// ks
	unsigned int sz_kr = ks_real.finalVals.size(), sz_ki = ks_imag.finalVals.size();
	sz_k = MAX(sz_kr, sz_ki);
	ks.resize(sz_k);
	for (unsigned int i = 0; i < sz_k; ++i)
		ks[i] = 0.0;
	if (sz_kr > 0)
		for (unsigned int i = 0; i < sz_kr; ++i)
			ks[i].real(ks_real[i]);
	if (sz_ki > 0)
		for (unsigned int i = 0; i < sz_ki; ++i)
			ks[i].imag(ks_imag[i]);
	// omegas
	unsigned int sz_omegar = omegas_real.finalVals.size(), sz_omegai = omegas_imag.finalVals.size();
	sz_omega = MAX(sz_omegar, sz_omegai);
	omegas.resize(sz_omega);
	for (unsigned int i = 0; i < sz_omega; ++i)
		omegas[i] = 0.0;
	if (sz_omegar > 0)
		for (unsigned int i = 0; i < sz_omegar; ++i)
			omegas[i].real(omegas_real[i]);
	if (sz_omegai > 0)
		for (unsigned int i = 0; i < sz_omegai; ++i)
			omegas[i].imag(omegas_imag[i]);
}

void Configuration::Compute_StaticKaFSystem()
{
	if (solutionMode != smt_MatsF_Static)
		return;
	double valLeft, valRight;
	Compute_left_right_prescribedVals(0.0, valLeft, valRight);
	globalF.resize(nfdof_domain);
	ZEROVEC(globalF, nfdof_domain);

	// Dirichlet, Neumann BCs
	if (isDG)
	{
		// Dirichlet BC
		for (unsigned int dbci = 0; dbci < num_DirichletBCs; ++dbci)
		{
			DCVECTOR interfaceF;
			interfaceF.resize(ndof_element);
			ZEROVEC(interfaceF, ndof_element);
			bool leftDomainInterface = (dbci == 0);
			double wBar = valRight;
			if (leftDomainInterface)
				wBar = valLeft;
			Compute_DirichletBoundaryInterface_ForceVector(*(e_DirichletBCs[dbci]),
				leftDomainInterface, wBar, interfaceF);

			unsigned int pos = eIndex_DirichletBCs[dbci];
			unsigned int st = element_start_dof[pos];

			for (unsigned int i = 0; i < ndof_element; ++i)
				globalF(i + st) += interfaceF(i);
		}
		// Neumann BC
		for (unsigned int nbci = 0; nbci < num_NeumannBCs; ++nbci)
		{
			DCVECTOR interfaceF;
			bool leftDomainInterface = (nbci == 0);
			double sigmaBarn = valRight;
			if (leftDomainInterface)
				sigmaBarn = valLeft;
			Compute_NeumannBoundaryInterface_ForceVector(*(e_NeumannBCs[nbci]),
				leftDomainInterface, sigmaBarn, interfaceF);

			unsigned int pos = eIndex_NeumannBCs[nbci];
			unsigned int st = element_start_dof[pos];

			for (unsigned int i = 0; i < ndof_element; ++i)
				globalF(i + st) += interfaceF(i);
		}
	}
	// if it has interior source term, we'll add the computation here

	// system solution
#if VCPP
	ZEROVEC(global_a, nfdof_domain);
	for (unsigned int i = 0; i < nfdof_domain; ++i)
		global_a(i) = globalF(i);
	DCMATRIX globalK_cpy;
	ZEROMAT(globalK_cpy, nfdof_domain, nfdof_domain);
	globalK_cpy = globalK;
	LUsolve(globalK_cpy, global_a);
#else
    // eigen-based Ka = F --- you may need to have 3 versions of DCVECTOR 1 VCPP, 2 otherwise

	// 1. Eigen-based Ka = F
	ZEROVEC(global_a, nfdof_domain);
	global_a = globalK.fullPivLu().solve(globalF);

#endif
}

void Configuration::Compute_BlochModeAnalysis()
{
	if (solutionMode != smt_Mats_OnlyBloch)
		return;

	string str_ser;
	toString(serialNumber, str_ser);
	string filename = outputFileNameWOExt + "_" + str_ser + "_Bloch.txt";
	fstream out(filename.c_str(), ios::out);
	out << setprecision(22);
	out << "sz_k\t" << sz_k << '\n';

	Dcomplex I(0.0, 1.0);

	for (unsigned ik = 0; ik < sz_k; ++ik)
	{
		//! 1. Forming k, gamma = exp(ikL), gamma_inv
		Dcomplex k = ks[ik];
		Dcomplex exp_v = domain_length * (I * k);
		gamma = exp(exp_v);
		gamma_inv = 1.0 / gamma;

		DCMATRIX MBloch, KBloch, CBloch;

		if (isDG)
		{ // DG start
			MBloch = globalM;
			KBloch = globalK;
			if (b_hasNonZeroDampingMatrix)
				CBloch = globalC;

			//! 2DG: Updating K for Bloch BC
			unsigned int interface_dof = 2 * ndof_element;
			unsigned interfacei = num_elements;
			const OneDimensionalElement *left_ePtr, *right_ePtr;
			int left_e_dof_st, right_e_dof_st;
			bool insideDomainInterface = false;
			left_e_dof_st = element_start_dof[interfacei - 1];
			left_ePtr = &elements[interfacei - 1];
			right_e_dof_st = element_start_dof[0];
			right_ePtr = &elements[0];
			DCMATRIX interfaceK, interfaceC;
			ZEROMAT(interfaceK, interface_dof, interface_dof);
			ZEROMAT(interfaceC, interface_dof, interface_dof);
			bool hasC = Compute_InteriorInterface_Matrices(*left_ePtr, *right_ePtr, insideDomainInterface,
				interfaceK, interfaceC);

			for (unsigned int i = 0; i < ndof_element; ++i)
			{
				for (unsigned int j = 0; j < ndof_element; ++j)
				{
					// RR
					KBloch(i + right_e_dof_st, j + right_e_dof_st) += interfaceK(i + ndof_element, j + ndof_element);
					// RL
					KBloch(i + right_e_dof_st, j + left_e_dof_st) += interfaceK(i + ndof_element, j);
					// LR
					KBloch(i + left_e_dof_st, j + right_e_dof_st) += interfaceK(i, j + ndof_element);
					// LL
					KBloch(i + left_e_dof_st, j + left_e_dof_st) += interfaceK(i, j);
				}
			}
			if (hasC)
			{
				for (unsigned int i = 0; i < ndof_element; ++i)
				{
					for (unsigned int j = 0; j < ndof_element; ++j)
					{
						// RR
						CBloch(i + right_e_dof_st, j + right_e_dof_st) += interfaceC(i + ndof_element, j + ndof_element);
						// RL
						CBloch(i + right_e_dof_st, j + left_e_dof_st) += interfaceC(i + ndof_element, j);
						// LR
						CBloch(i + left_e_dof_st, j + right_e_dof_st) += interfaceC(i, j + ndof_element);
						// LL
						CBloch(i + left_e_dof_st, j + left_e_dof_st) += interfaceC(i, j);
					}
				}
			}
		} // DG end
		else // CFEM start
		{
			// For Ali
			// use the relation that last dof = gamma first dof to shrink global M, and K by size 1 to their Bloch form
		} // CFEM end

		//!. 3 Bloch Calculations - Now the matrices are simply output

		out << "ik\t" << ik << '\t';
		out << "k\t" << k << '\t';
		out << "gamma\t" << gamma << '\t';
		out << "gamma_inv\t" << gamma_inv << '\n';
		out << "KBloch\n" << KBloch << '\n';
		out << "MBloch\n" << MBloch << '\n';
		if (b_hasNonZeroDampingMatrix)
			out << "CBloch\n" << CBloch << '\n';

		// For Ali
		//! 4. You can do Bloch calculations here or ideally read this matrices outside and do the calculation there. 
		// In the latter case it may be a good idea to write the matrices in binary
	}
}

void Configuration::Process_Output_GlobalMatrices()
{
	string str_ser;
	toString(serialNumber, str_ser);
    string filename = outputFileNameWOExt + "_" + str_ser + ".txt";
    string filename_k = outputFileNameWOExt + "_k_matrix" + str_ser + ".txt";
    string filename_m = outputFileNameWOExt + "_m_matrix" + str_ser + ".txt";
    string filename_c = outputFileNameWOExt + "_c_matrix" + str_ser + ".txt";
    string filename_eigen = outputFileNameWOExt + "_eigen_" + str_ser + ".txt";
    fstream out(filename.c_str(), ios::out);
    fstream out_k(filename_k.c_str(), ios::out);
    fstream out_m(filename_m.c_str(), ios::out);
    fstream out_c(filename_c.c_str(), ios::out);
    fstream out_eigen(filename_eigen.c_str(), ios::out);
    out << "wf_type\t" << wf_type << '\n';
    out << "polyOrder\t" << polyOrder << '\n';
    out << "ndof_domain\t" << ndof_domain << '\n';
    out << "star option\t" << sOption << '\n';
    out << "dg_eps\t" << dg_eps << '\n';
    out << "etaI\t" << etaI << '\n';
    out << "etaB\t" << etaB << '\n';
    //out << "K\n" << globalK << '\n';
    out << "K.sum\t" << globalK.sum() << '\n';
    //out << "M\n" << globalM << '\n';
    out << "M.sum\t" << globalM.sum() << '\n';
    out << "b_hasNonZeroDampingMatrix\t" << b_hasNonZeroDampingMatrix << '\n';
    //out << "C\n" << globalC << '\n';
    out << "C.sum\n" << globalC.sum() << '\n';
    Eigen::IOFormat OctaveFmt(Eigen::FullPrecision, 0, ", ", " \n", "", "", " ", " ");
    out_k << globalK.format(OctaveFmt) << '\n';
    out_m << globalM.format(OctaveFmt) << '\n';
    out_c << globalC.format(OctaveFmt) << '\n';
	if (solutionMode == smt_MatsF_Static)
	{
		out << "globalF\n" << globalF << '\n';
		out << "global_a\n" << global_a << '\n';
	}
	if (solutionMode == smt_Mats_OnlyNaturalMode)
	{
#if !VCPP       
#if USE_COMPLEX
		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
        ces.compute (globalM.inverse()*globalK, /* computeEigenvectors = */ false);
        // print stuff
        out << "eigenvalues\n" << ces.eigenvalues() << '\n';
#else
        Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;
        ces.compute (globalM.inverse()*globalK, /* computeEigenvectors = */ false);
        // print stuff
        out << "eigenvalues\n" << ces.eigenvalues() << '\n';
        out_eigen << ces.eigenvalues().format(OctaveFmt) << '\n';
        if (b_hasNonZeroDampingMatrix)
        {
            Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces2;

            Eigen::MatrixXd M_bar(2*nfdof_domain,2*nfdof_domain), K_bar(2*nfdof_domain,2*nfdof_domain);
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nfdof_domain,nfdof_domain);
            Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(nfdof_domain,nfdof_domain);
            M_bar << I , Z, Z , globalM;
            K_bar << Z , -I, globalK, globalC;
            ces2.compute (M_bar.inverse()*K_bar, /* computeEigenvectors = */ false);
            out << "eigenvalues2\n" << ces2.eigenvalues() << '\n';
            out_eigen << ces2.eigenvalues().format(OctaveFmt) << '\n';
        }

#endif
#endif
	//
		//out << "eigenvalues\n" << ges.eigenvalues() << '\n';
	//
		//Eigen::EigenSolver<Eigen::MatrixXd> es;
	//
		//es.compute (globalM.inverse()*globalK, /* computeEigenvectors = */ false);
	//
		//out << "eigenvalues\n" << es.eigenvalues() << '\n';
	}
#if USE_COMPLEX
#if !VCPP
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;

    ces.compute (globalM.inverse()*globalK);

    out << "eigenvalues\n" << ces.eigenvalues() << '\n';
#endif
#endif
}


StarW1s::StarW1s()
{
    ss_f_sigmaL = 0.0;
    ss_f_sigmaR = 0.0;
    ss_f_wL = 0.0;
    ss_f_wR = 0.0;
    ws_f_sigmaL = 0.0;
    ws_f_sigmaR = 0.0;
    ws_f_wL = 0.0;
    ws_f_wR = 0.0;
  
}

double ComputeLoad(double t, int loadNumber, const vector<double>& loadValues)
{
	if (loadNumber == 0)
		return loadValues[0];
	if (loadNumber == 1)
		return loadValues[0] + loadValues[1] * t;
	return 0.0;
}
