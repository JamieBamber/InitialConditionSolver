/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SetLevelData.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "MyPhiFunction.H"
#include "MyPotentialFunction.H"
#include "PoissonParameters.H"
#include "SetAij.H"
#include "SetLevelDataCpp.H"
//#include "SetLevelDataF_F.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij, scalar field phi
// and \bar Aij = psi^2 A_ij.
// For now the default setup is 2 Bowen York BHs plus a scalar field
// with some initial user specified configuration
void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars,
                            LevelData<FArrayBox> &a_dpsi, const RealVect &a_dx,
                            const PoissonParameters &a_params) {

  CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

  DataIterator dit = a_multigrid_vars.dataIterator();
  const DisjointBoxLayout &grids = a_multigrid_vars.disjointBoxLayout();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
    FArrayBox &dpsi_box = a_dpsi[dit()];
    Box b = multigrid_vars_box.box();
    Box b_no_ghosts = grids[dit()];
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit) {

      // work out location on the grid
      IntVect iv = bit();
      RealVect loc(iv + 0.5 * RealVect::Unit);
      loc *= a_dx;
      loc -= a_params.domainLength / 2.0;

      // note that we don't include the singular part of psi
      // for the BHs - this is added at the output data stage
      // and when we calculate psi_0 in the rhs etc
      // as it already satisfies Laplacian(psi) = 0
      multigrid_vars_box(iv, c_psi_0) = a_params.psi_0;
      multigrid_vars_box(iv, c_V0_0) = 0.0;
      multigrid_vars_box(iv, c_V1_0) = 0.0;
      multigrid_vars_box(iv, c_V2_0) = 0.0;

      dpsi_box(iv, c_psi) = 0.0;
      dpsi_box(iv, c_V0) = 0.0;
      dpsi_box(iv, c_V1) = 0.0;
      dpsi_box(iv, c_V2) = 0.0;

      // set phi and pi according to user defined function
      multigrid_vars_box(iv, c_phi_0) = my_phi_function(loc, a_params);

      multigrid_vars_box(iv, c_pi_0) = my_pi_function(loc, a_params);
    }
  }
} // end set_initial_conditions

// set the rhs source for the poisson eqn
void set_rhs(LevelData<FArrayBox> &a_rhs,
             LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx,
             const PoissonParameters &a_params) {

  CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

  DataIterator dit = a_rhs.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
    FArrayBox &rhs_box = a_rhs[dit()];
    rhs_box.setVal(0.0, 0);
    Box this_box = rhs_box.box(); // no ghost cells

    // calculate the laplacian of Psi across the box
    FArrayBox laplace_multigrid(this_box, NUM_CONSTRAINTS_VARS);
    get_laplacian(this_box, multigrid_vars_box, Interval(c_psi, c_V2), a_dx,
                  laplace_multigrid, a_params);

    // calculate gradients for constructing rho and Aij
    FArrayBox grad_multigrid(this_box, 3 * NUM_MULTIGRID_VARS);
    get_grad(this_box, multigrid_vars_box, Interval(c_psi_0, c_pi_0), a_dx,
             grad_multigrid, a_params);

    FArrayBox grad2_multigrid(this_box, 3 * NUM_MULTIGRID_VARS);
    get_grad2(this_box, multigrid_vars_box, Interval(c_psi_0, c_pi_0), a_dx,
              grad2_multigrid, a_params);

    FArrayBox mixed_grad2_multigrid(this_box, 3 * NUM_MULTIGRID_VARS);
    get_mixed_grad2(this_box, multigrid_vars_box, Interval(c_psi_0, c_pi_0),
                    a_dx, mixed_grad2_multigrid, a_params);

    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect loc(iv + 0.5 * RealVect::Unit);
      loc *= a_dx;
      loc -= a_params.domainLength / 2.0;

      // rhs = term/8 psi_0^5 - 2 pi rho_grad psi_0  - laplacian(psi_0)
      // JCAurre: Calculate potential and derivative
      Real V_of_phi, dVdphi;
      my_potential_function(V_of_phi, dVdphi, multigrid_vars_box(iv, c_phi_0),
                            a_params);

      // JCAurre: useful
      Real psi_bh = set_binary_bh_psi(loc, a_params);
      Real psi_0 = multigrid_vars_box(iv, c_psi_0) + psi_bh;
      Real K_0 = multigrid_vars_box(iv, c_K_0);
      Real pi_0 = multigrid_vars_box(iv, c_pi_0);

      Real d_phi[3], dd_phi[3][3], d_pi[3], Aij[3][3], d_Aij[3][3][3];
      set_deriv1(d_phi, iv, grad_multigrid, c_phi_0);
      set_deriv2(dd_phi, iv, grad2_multigrid, mixed_grad2_multigrid, c_phi_0);
      set_deriv1(d_pi, iv, grad_multigrid, c_pi_0);

      set_Aij_0(Aij, multigrid_vars_box, iv, loc, a_dx, a_params,
                grad_multigrid);

      set_deriv_Aij_0(d_Aij, multigrid_vars_box, iv, loc, a_dx, a_params,
                      grad2_multigrid, mixed_grad2_multigrid);

      set_binary_bh_Aij(multigrid_vars_box, iv, loc, a_params);

      // JCAurre: compute rhograd from gradients of phi
      Real rho_gradient = 0;
      for (int i = 0; i < SpaceDim; i++) {
        rho_gradient += 0.5 * d_phi[i] * d_phi[i];
      }

      // Also \bar  A_ij \bar A^ij
      Real A2 = 0.0;
      for (int i = 0; i < SpaceDim; i++) {
        for (int j = 0; j < SpaceDim; j++) {
          A2 += Aij[i][j] * Aij[i][j];
        }
      }

      Real Aij_d_Aij[3] = {0, 0, 0};
      for (int i = 0; i < SpaceDim; i++) {
        for (int j = 0; j < SpaceDim; j++) {
          for (int k = 0; k < SpaceDim; k++) {
            Aij_d_Aij[i] += Aij[j][k] * d_Aij[j][k][i];
          }
        }
      }

      rhs_box(iv, c_psi) =
          0.125 *
              (2.0 / 3.0 * K_0 * K_0 - 8.0 * M_PI * a_params.G_Newton *
                                           (pow(pi_0, 2.0) + 2.0 * V_of_phi)) *
              pow(psi_0, 5.0) -
          0.125 * A2 * pow(psi_0, -7.0) -
          2.0 * M_PI * a_params.G_Newton * rho_gradient * psi_0 -
          laplace_multigrid(iv, c_psi);

      // JCAurre: Added rhs for new constraint variables.
      rhs_box(iv, c_V0) = -8.0 * M_PI * pow(psi_0, 6.0) * pi_0 * d_phi[0] -
                          laplace_multigrid(iv, c_V0);
      rhs_box(iv, c_V1) = -8.0 * M_PI * pow(psi_0, 6.0) * pi_0 * d_phi[1] -
                          laplace_multigrid(iv, c_V1);
      rhs_box(iv, c_V2) = -8.0 * M_PI * pow(psi_0, 6.0) * pi_0 * d_phi[2] -
                          laplace_multigrid(iv, c_V2);

      // JCAurre: Add grad K terms
      Real gradphi_grad2phi[3] = {0, 0, 0};
      for (int i = 0; i < SpaceDim; i++) {
        for (int j = 0; j < SpaceDim; j++) {
          gradphi_grad2phi[i] += d_phi[j] * dd_phi[j][i];
        }
      }

      Real min_K_0 = -1.e-15;
      if (abs(K_0) < abs(min_K_0)) {
        K_0 = min_K_0;
      }

      RealVect gradK;
      for (int i = 0; i < SpaceDim; i++) {
        gradK[i] =
            pow(K_0, -1.0) *
            (12.0 * M_PI *
                 (pi_0 * d_pi[i] + pow(psi_0, -4.0) * gradphi_grad2phi[i] +
                  dVdphi * d_phi[i]) +
             1.5 * pow(psi_0, -12.0) * Aij_d_Aij[i]);
      }

      rhs_box(iv, c_V0) += 2. / 3. * pow(psi_0, 6.0) * gradK[0];
      rhs_box(iv, c_V1) += 2. / 3. * pow(psi_0, 6.0) * gradK[1];
      rhs_box(iv, c_V2) += 2. / 3. * pow(psi_0, 6.0) * gradK[2];
    }
  }
} // end set_rhs

// Set the integrand for the integrability condition for constant K
// when periodic BCs set
void set_integrability(LevelData<FArrayBox> &a_integrand,
                       LevelData<FArrayBox> &a_multigrid_vars,
                       const RealVect &a_dx,
                       const PoissonParameters &a_params) {

  CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

  DataIterator dit = a_integrand.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
    FArrayBox &integrand_box = a_integrand[dit()];
    integrand_box.setVal(0.0, 0);
    Box this_box = integrand_box.box(); // no ghost cells

    // calculate the laplacian of Psi across the box
    FArrayBox laplace_multigrid(this_box, NUM_CONSTRAINTS_VARS);
    get_laplacian(this_box, multigrid_vars_box, Interval(c_psi, c_V2), a_dx,
                  laplace_multigrid, a_params);

    // calculate gradients for constructing rho and Aij
    FArrayBox grad_multigrid(this_box, 3 * NUM_MULTIGRID_VARS);
    get_grad(this_box, multigrid_vars_box, Interval(c_psi_0, c_phi_0), a_dx,
             grad_multigrid, a_params);

    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect loc(iv + 0.5 * RealVect::Unit);
      loc *= a_dx;
      loc -= a_params.domainLength / 2.0;

      // integrand = -1.5*term + 1.5 * \bar A_ij \bar A^ij psi_0^-12 +
      // 24 pi rho_grad psi_0^-4  + 12*laplacian(psi_0)*psi^-5
      // JCAurre: Calculate potential and derivative
      Real V_of_phi, dVdphi;
      my_potential_function(V_of_phi, dVdphi, multigrid_vars_box(iv, c_phi_0),
                            a_params);

      // JCAurre: useful
      Real psi_bh = set_binary_bh_psi(loc, a_params);
      Real psi_0 = multigrid_vars_box(iv, c_psi_0) + psi_bh;
      Real pi_0 = multigrid_vars_box(iv, c_pi_0);

      Real d_phi[3];
      set_deriv1(d_phi, iv, grad_multigrid, c_phi_0);
      Real Aij[3][3];
      set_Aij_0(Aij, multigrid_vars_box, iv, loc, a_dx, a_params,
                grad_multigrid);
      set_binary_bh_Aij(multigrid_vars_box, iv, loc, a_params);

      // JCAurre: compute rhograd from gradients of phi
      Real rho_gradient = 0;
      for (int i = 0; i < SpaceDim; i++) {
        rho_gradient += 0.5 * d_phi[i] * d_phi[i];
      }

      // Also \bar  A_ij \bar A^ij
      Real A2 = 0.0;
      for (int i = 0; i < SpaceDim; i++) {
        for (int j = 0; j < SpaceDim; j++) {
          A2 += Aij[i][j] * Aij[i][j];
        }
      }

      Real K_0_sq =
          -1.5 * (0 - 8.0 * M_PI * a_params.G_Newton *
                          (pow(pi_0, 2.0) + 2.0 * V_of_phi)) +
          1.5 * A2 * pow(psi_0, -12.0) +
          24.0 * M_PI * a_params.G_Newton * rho_gradient * pow(psi_0, -4.0) +
          12.0 * laplace_multigrid(iv, c_psi) * pow(psi_0, -5.0);

      integrand_box(iv, c_psi) =
          -1.5 * (2.0 / 3.0 * K_0_sq - 8.0 * M_PI * a_params.G_Newton *
                                           (pow(pi_0, 2.0) + 2.0 * V_of_phi)) +
          1.5 * A2 * pow(psi_0, -12.0) +
          24.0 * M_PI * a_params.G_Newton * rho_gradient * pow(psi_0, -4.0) +
          12.0 * laplace_multigrid(iv, c_psi) * pow(psi_0, -5.0);
      integrand_box(iv, c_V0) =
          -8.0 * M_PI * pow(psi_0, 6.0) * pi_0 * d_phi[0] -
          laplace_multigrid(iv, c_V0);
      integrand_box(iv, c_V1) =
          -8.0 * M_PI * pow(psi_0, 6.0) * pi_0 * d_phi[1] -
          laplace_multigrid(iv, c_V1);
      integrand_box(iv, c_V2) =
          -8.0 * M_PI * pow(psi_0, 6.0) * pi_0 * d_phi[2] -
          laplace_multigrid(iv, c_V2);

      // Set value for K
      multigrid_vars_box(iv, c_K_0) =
          -sqrt(K_0_sq); // be careful when K=0, maybe discontinuity
    }
  }
} // end set_constant_K_integrand

// set the regrid condition - abs value of this drives AMR regrid
void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                          LevelData<FArrayBox> &a_multigrid_vars,
                          const RealVect &a_dx,
                          const PoissonParameters &a_params) {

  CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

  DataIterator dit = a_condition.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
    FArrayBox &condition_box = a_condition[dit()];
    condition_box.setVal(0.0, 0);
    Box this_box = condition_box.box(); // no ghost cells

    // calculate gradients for constructing rho and Aij
    FArrayBox grad_multigrid(this_box, 3 * NUM_MULTIGRID_VARS);
    get_grad(this_box, multigrid_vars_box, Interval(c_psi_0, c_phi_0), a_dx,
             grad_multigrid, a_params);

    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect loc(iv + 0.5 * RealVect::Unit);
      loc *= a_dx;
      loc -= a_params.domainLength / 2.0;

      // JCAurre: Calculate potential and derivative
      Real V_of_phi, dVdphi;
      my_potential_function(V_of_phi, dVdphi, multigrid_vars_box(iv, c_phi_0),
                            a_params);

      // JCAurre: useful
      Real psi_bh = set_binary_bh_psi(loc, a_params);
      Real psi_0 = multigrid_vars_box(iv, c_psi_0) + psi_bh;
      Real K_0 = multigrid_vars_box(iv, c_K_0);
      Real pi_0 = multigrid_vars_box(iv, c_pi_0);

      Real d_phi[3];
      set_deriv1(d_phi, iv, grad_multigrid, c_phi_0);
      Real Aij[3][3];
      set_Aij_0(Aij, multigrid_vars_box, iv, loc, a_dx, a_params,
                grad_multigrid);
      set_binary_bh_Aij(multigrid_vars_box, iv, loc, a_params);

      // JCAurre: compute rhograd from gradients of phi
      Real rho_gradient = 0;
      for (int i = 0; i < SpaceDim; i++) {
        rho_gradient += 0.5 * d_phi[i] * d_phi[i];
      }
      // Also \bar  A_ij \bar A^ij
      Real A2 = 0.0;
      for (int i = 0; i < SpaceDim; i++) {
        for (int j = 0; j < SpaceDim; j++) {
          A2 += Aij[i][j] * Aij[i][j];
        }
      }

      // the condition is similar to the rhs but we take abs
      // value of the contributions and add in BH criteria
      condition_box(iv, 0) =
          1.5 * abs((2.0 / 3.0 * K_0 * K_0 +
                     8.0 * M_PI * a_params.G_Newton *
                         (pow(pi_0, 2.0) + 2.0 * V_of_phi))) +
          1.5 * A2 * pow(psi_0, -7.0) +
          24.0 * M_PI * a_params.G_Newton * abs(rho_gradient) *
              pow(psi_0, 1.0) +
          log(psi_0);

      // JCAurre: Maybe one wants to add regridding in the momentum
      // constraints
     // condition_box(iv, c_V0) = 0.;
     // condition_box(iv, c_V1) = 0.;
     // condition_box(iv, c_V2) = 0.;
    }
  }
} // end set_regrid_condition

// Add the correction to psi0 after the solver operates
void set_update_psi0(LevelData<FArrayBox> &a_multigrid_vars,
                     LevelData<FArrayBox> &a_dpsi,
                     const Copier &a_exchange_copier) {

  // first exchange ghost cells for dpsi so they are filled with the correct
  // values
  a_dpsi.exchange(a_dpsi.interval(), a_exchange_copier);

  DataIterator dit = a_multigrid_vars.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
    FArrayBox &dpsi_box = a_dpsi[dit()];

    Box this_box = multigrid_vars_box.box();
    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();

      // JCAurre: update constraint variables
      multigrid_vars_box(iv, c_psi_0) += dpsi_box(iv, c_psi);
      multigrid_vars_box(iv, c_V0_0) += dpsi_box(iv, c_V0);
      multigrid_vars_box(iv, c_V1_0) += dpsi_box(iv, c_V1);
      multigrid_vars_box(iv, c_V2_0) += dpsi_box(iv, c_V2);
    }
  }
}

// The coefficient of the I operator on dpsi
void set_a_coef(LevelData<FArrayBox> &a_aCoef,
                LevelData<FArrayBox> &a_multigrid_vars,
                const PoissonParameters &a_params, const RealVect &a_dx) {

  CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

  DataIterator dit = a_aCoef.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &aCoef_box = a_aCoef[dit()];

    // JCAurre: aCoef is now set to zero as Ham is algebraic
    // equation and Mom is linear
    for (int iconstraint = 0; iconstraint < NUM_CONSTRAINTS_VARS;
         iconstraint++) {
      aCoef_box.setVal(0.0, iconstraint);
    }
  }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
// the rhs source of the Poisson eqn
void set_b_coef(LevelData<FArrayBox> &a_bCoef,
                const PoissonParameters &a_params, const RealVect &a_dx) {

  CH_assert(a_bCoef.nComp() == NUM_CONSTRAINTS_VARS);
  int comp_number = 0;

  for (DataIterator dit = a_bCoef.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &bCoef_box = a_bCoef[dit()];
    // JCAurre: Loop to set bCoef=1 for all constraint variables
    for (int iconstraint = 0; iconstraint < NUM_CONSTRAINTS_VARS;
         iconstraint++) {
      bCoef_box.setVal(1.0, iconstraint);
    }
  }
}

// used to set output data for all ADM Vars for GRChombo restart
void set_output_data(LevelData<FArrayBox> &a_grchombo_vars,
                     LevelData<FArrayBox> &a_multigrid_vars,
                     const PoissonParameters &a_params, const RealVect &a_dx) {

  CH_assert(a_grchombo_vars.nComp() == NUM_GRCHOMBO_VARS);
  CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

  const DisjointBoxLayout &grids = a_grchombo_vars.disjointBoxLayout();
  DataIterator dit = a_grchombo_vars.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &grchombo_vars_box = a_grchombo_vars[dit()];
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

    // first set everything to zero
    for (int comp = 0; comp < NUM_GRCHOMBO_VARS; comp++) {
      grchombo_vars_box.setVal(0.0, comp);
    }

    // now set non zero terms - const across whole box
    // Conformally flat, and lapse = 1
    grchombo_vars_box.setVal(1.0, c_h11);
    grchombo_vars_box.setVal(1.0, c_h22);
    grchombo_vars_box.setVal(1.0, c_h33);
    grchombo_vars_box.setVal(1.0, c_lapse);

    // now non constant terms by location
    Box this_box = grchombo_vars_box.box();
    BoxIterator bit(this_box);
    Box this_box_ng = grids[dit()];

    // calculate gradients for constructing Aij
    FArrayBox grad_multigrid(this_box_ng, 3 * NUM_MULTIGRID_VARS);
    get_grad(this_box_ng, multigrid_vars_box, Interval(c_V0_0, c_V2_0), a_dx,
             grad_multigrid, a_params);

    // Aij is defined to be zero in the ghost cells, so be careful with
    // fixed BCs
    BoxIterator bit_no_ghosts(this_box_ng);
    for (bit_no_ghosts.begin(); bit_no_ghosts.ok(); ++bit_no_ghosts) {

      // work out location on the grid
      IntVect iv = bit_no_ghosts();

      // set the phi value - need the distance from centre
      RealVect loc(iv + 0.5 * RealVect::Unit);
      loc *= a_dx;
      loc -= a_params.domainLength / 2.0;

      // JCAurre: set Aij from components of vector W
      Real Aij[3][3]; // not used afterwards
      set_Aij_0(Aij, multigrid_vars_box, iv, loc, a_dx, a_params,
                grad_multigrid);
    }

    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect loc(iv + 0.5 * RealVect::Unit);
      loc *= a_dx;
      loc -= a_params.domainLength / 2.0;

      // GRChombo conformal factor chi = psi^-4
      Real psi_bh = set_binary_bh_psi(loc, a_params);
      Real chi = pow(multigrid_vars_box(iv, c_psi) + psi_bh, -4.0);
      grchombo_vars_box(iv, c_chi) = chi;
      Real factor = pow(chi, 1.5);

      // Copy phi and Aij across - note this is now \tilde Aij not \bar
      // Aij
      grchombo_vars_box(iv, c_phi) = multigrid_vars_box(iv, c_phi_0);
      grchombo_vars_box(iv, c_Pi) = multigrid_vars_box(iv, c_pi_0);

      grchombo_vars_box(iv, c_K) = multigrid_vars_box(iv, c_K_0);
      grchombo_vars_box(iv, c_A11) = multigrid_vars_box(iv, c_A11_0) * factor;
      grchombo_vars_box(iv, c_A12) = multigrid_vars_box(iv, c_A12_0) * factor;
      grchombo_vars_box(iv, c_A13) = multigrid_vars_box(iv, c_A13_0) * factor;
      grchombo_vars_box(iv, c_A22) = multigrid_vars_box(iv, c_A22_0) * factor;
      grchombo_vars_box(iv, c_A23) = multigrid_vars_box(iv, c_A23_0) * factor;
      grchombo_vars_box(iv, c_A33) = multigrid_vars_box(iv, c_A33_0) * factor;
    }
  }
}
