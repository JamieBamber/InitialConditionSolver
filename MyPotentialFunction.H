/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYPOTENTIALFUNCTION_HPP
#define MYPOTENTIALFUNCTION_HPP

void my_potential_function(Real &V_of_phi, Real &dVdphi_Re, Real &dVdphi_Im, const Real &phi_Re_here, const Real &phi_Im_here,
                           const PoissonParameters &a_params)
{
    // Minimally coupled scalar field
    V_of_phi = 0.5 * (pow(a_params.pot_mu * phi_Re_here, 2) + pow(a_params.pot_mu * phi_Im_here, 2));
    dVdphi_Re = pow(a_params.pot_mu,2) * phi_Re_here;
    dVdphi_Im = pow(a_params.pot_mu,2) * phi_Im_here;
}

#endif /* MYPOTENTIALFUNCTION_HPP */
