/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYMATTERFUNCTIONS_HPP
#define MYMATTERFUNCTIONS_HPP

#include "SetPsiAndAij.H"

Real my_phi_Re_function(const RealVect loc, const PoissonParameters &a_params) {

    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return a_params.phi_amplitude *
             exp(-pow(rr * a_params.phi_kappa, 2.0));
}

Real my_phi_Im_function(const RealVect loc, const PoissonParameters &a_params) {

   Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
   return 0.0;
}

Real my_Pi_Re_function(const RealVect loc, const PoissonParameters &a_params) {

    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return 0.0;
}

Real my_Pi_Im_function(const RealVect loc, const PoissonParameters &a_params) {

    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return - a_params.scalar_omega * a_params.phi_amplitude *
             exp(-pow(rr * a_params.phi_kappa, 2.0));
}

Real my_potential_function(const Real &phi_Re_here, const Real &phi_Im_here,
                           const PoissonParameters &a_params)
{
    return a_params.mu * a_params.mu * (phi_Re_here * phi_Re_here + phi_Im_here * phi_Im_here);
}

Real my_rho_function(const IntVect a_iv, const RealVect &a_dx,
                     FArrayBox &a_multigrid_vars_box,
                     const PoissonParameters &a_params)
{

    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh = set_binary_bh_psi(loc, a_params);
    Real psi_0 = psi_reg + psi_bh;

    Real V_of_phi =
        my_potential_function(a_multigrid_vars_box(a_iv, c_phi_Re_0), a_multigrid_vars_box(a_iv, c_phi_Im_0), a_params);
    Real Pi_Re_0 = a_multigrid_vars_box(a_iv, c_Pi_Re_0);
    Real Pi_Im_0 = a_multigrid_vars_box(a_iv, c_Pi_Im_0);
    Tensor<1, Real, SpaceDim> d1_phi_Re =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_Re_0);
    Tensor<1, Real, SpaceDim> d1_phi_Im =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_Im_0);
    Real d1_phi_squared = 0;
    FOR1(i) { d1_phi_squared += d1_phi_Re[i] * d1_phi_Re[i] + d1_phi_Im[i] * d1_phi_Im[i]; }
    return 0.5 * (Pi_Re_0 * Pi_Re_0 + Pi_Im_0 * Pi_Im_0) + V_of_phi +
           0.5 * d1_phi_squared * pow(psi_0, -4.0);
}

Real my_Si_function(Tensor<1, Real> &Si, const IntVect a_iv,
                    const RealVect &a_dx, FArrayBox &a_multigrid_vars_box,
                    const PoissonParameters &a_params)
{
    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    Real Pi_Re_0 = a_multigrid_vars_box(a_iv, c_Pi_Re_0);
    Real Pi_Im_0 = a_multigrid_vars_box(a_iv, c_Pi_Im_0);
    Tensor<1, Real, SpaceDim> d1_phi_Re =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_Re_0);
    Tensor<1, Real, SpaceDim> d1_phi_Im =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_Im_0);

    FOR1(i) { Si[i] = -Pi_Re_0 * d1_phi_Re[i] - Pi_Im_0 * d1_phi_Im[i]; }
}

#endif /* MYMATTERFUNCTIONS_HPP */
