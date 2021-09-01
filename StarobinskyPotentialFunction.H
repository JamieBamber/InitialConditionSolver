/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYPOTENTIALFUNCTION_HPP
#define MYPOTENTIALFUNCTION_HPP

void my_potential_function(Real &V_of_phi, Real &dVdphi, const Real &phi_here,
                           const PoissonParameters &a_params)
{
    // Starobinsky
    Real mpl = 1.0;
    V_of_phi = a_params.pot_Lambda * pow(1. - exp((phi_here / a_params.pot_mu * pow(mpl, -1))), 2);
    dVdphi = -2. * a_params.pot_Lambda / a_params.pot_mu * (1 - exp((phi_here / a_params.pot_mu * pow(mpl, -1)))) *
             exp((phi_here / a_params.pot_mu * pow(mpl, -1))) * pow(mpl, -1);
}

#endif /* MYPOTENTIALFUNCTION_HPP */
