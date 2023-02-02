#ifndef CONFIGURATION__H
#define CONFIGURATION__H

#include "util/OneDimensionalElement.h"
#include "util/globalMacros.h"
#include "util/globalTypesClasses.h"

// weights of the star values for one side of the interface
// f_ -> factor of
// w: stands for v or u depending on the formulation
class StarW1s
{
  public:
    // sigma_nStar = ss_f_sigmaL * sigmaL + ss_f_sigmaR * sigmaR + ss_f_wL * wL +
    // ss_f_wR * wR
    NUMBR ss_f_sigmaL, ss_f_sigmaR, ss_f_wL, ss_f_wR;
    // w_Star = ws_f_sigmaL * sigmaL + ws_f_sigmaR * sigmaR + ws_f_wL * wL +
    // ws_f_wR * wR
    NUMBR ws_f_sigmaL, ws_f_sigmaR, ws_f_wL, ws_f_wR;
};

class StarW2s
{
  public:
    // side_weights[SDL]: the weights of the star values for the left side
    // side_weights[SDR]: // right side
    StarW1s side_weights[NUM_SIDES];
};

class Configuration
{
  public:
    Configuration();
    void Main_SolveDomain(string configName);

    // inputs
    unsigned int polyOrder;
    // whether we want the mass-lumpting option or not
    bool lumpMass;

    WeakFormulationType wf_type;
    // star options for DG method
    StarOption sOption;
    DG_epsilonT dg_eps;
    // for DG, non-Riemann star options sigma*n has a penalty multiplying (wR -
    // wL), that weight is computed from eta below etaI: for interior interfaces
    // etaB: for boundary faces
    double etaI, etaB;

    // computed from inputs
    bool isDG;
    unsigned int num_fields;

    // vHat (rho vDot + damping v - div. S - source) -> weight_BLM_is_velocity =
    // true uHat // -> weight_BLM_is_velocity = false
    bool weight_BLM_is_velocity;
    // what is the position of the weight relative to unknown fields of the
    // element
    unsigned int field_pos_weight_BLM;
    BoundaryConditionT leftBC, rightBC;
    bool isBlochModeAnalysis;
    Dcomplex gamma, gamma_inv; // gamma = exp(i k L) // k wavenumber, L is domain
                               // length -> needed for Bloch analysis

    bool isHelmholtz;
    double omega;

    /// large class members
    // we read / have already initialized [he, E, rho, damping] of elements[i].elementProps for each element and
    // calculate Z, c, etc. from this
    vector<OneDimensionalElement> elements; // for uniform case; elements.resize(num_elements)
    unsigned int num_elements;
    double xm, xM, E_uniform, rho_uniform, damping_uniform; // hE = (xM - xm)/num_elements
    ////////////////////////////////////////////
    // Computed
    OneDimensionalParentElement parentElement;
    unsigned int ndof_parent_element;
    unsigned int ndof_element;
    // f: in this context: all dofs that go to global matrices
    // for CFEM:
    //			f: free dofs (inside the domain, Neumann BC, ...)
    //			p: prescribed dofs (only Dirichlet BCs on the sides)
    // for DG: ndof_domain = nfdof_domain, npdof_domain = 0,
    //			"f" can correspond to Dirichlet and Neumann BCs as they go to global matrix
    unsigned int ndof_domain, nfdof_domain, npdof_domain;
    unsigned int number_of_nodes;
    // this is the entire domain nodal dof map: size = number_of_nodes
    vector<int> nodedof_map;
    // these are element dof maps
    vector<vector<int>> edof_2_globalDofMap;
    // these are indexed by ei (element index) and store the start and end dof of element (element goes grom  = start; <
    // end)
    vector<int> element_start_dof;
    vector<int> element_end_dof;

    /////////////////////////////////// Global, K, M, C
    // stores whether we need a damping matrix
    bool b_hasNonZeroDampingMatrix;
    DCMATRIX globalK, globalC, globalM;

  private:
    // functions called in Main_SolveDomain
    void Read_ConfigFile(string configName);
    void Initialize_Configuration();
    void Read_ElementGeometryProperties();
    void Form_ElementMatrices();
    void AssembleGlobalMatrices_DG(bool assembleMassIn);
    void AssembleGlobalMatrices_CFEM(bool assembleMassIn);
    // I think everything can be solved here, except Bloch analysis that should be written to a file for Ali
    void Process_Output_GlobalMatrices();

    /////////////////////////////////////////////////////

    inline bool IsWVelocity() const
    {
        return ((wf_type == DG_2FUV) || (wf_type == DG_1F_vStar));
    };

    // for 1F-formulation we end up with Ma'' + Ca' + Ka = F -> damping only if we have nonzero damping at element level
    // or we have 1Dvs formulation for 2F-formulation Ma' + Ka = F -> no damping If we have transmitting BC we need to
    // modify this function
    bool HasNonZeroDampingMatrix() const;

    ///	Auxiliary functions

    // inputs: left_ep, right_ep -> const makes it clearer that these are inputs
    //		   bool insideDomainInterface: periodic and Bloch BCs are
    // interior
    // interfaces. Bloch BC involves the value gamma 				other
    // interfaces that fall inside the domain are labeled as
    // "insideDomainInterface = true" for Bloch and periodic BC
    // insideDomainInterface = false outputs: twoSideWeights
    void Compute_DG_Star_Weights_4_Inteior_Interface(const ElementProperties &left_ep,
                                                     const ElementProperties &right_ep, bool insideDomainInterface,
                                                     StarW2s &twoSideWeights) const; // const here says that this object

    // (configuration) is not changing itself
    // return true if it has nonzero C
    bool Assemble_InteriorInterface_Matrices_2_Global_System(const OneDimensionalElement &left_e,
                                                             const OneDimensionalElement &right_e,
                                                             bool insideDomainInterface, DCMATRIX &interfaceK,
                                                             DCMATRIX &interfaceC) const;
};

#endif