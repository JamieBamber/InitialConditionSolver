/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYPOTENTIALFUNCTION_HPP
#define MYPOTENTIALFUNCTION_HPP

void my_potential_function(Real &V_of_phi, Real &dVdphi, const Real &phi_here,
                           const PoissonParameters &a_params)
{
    V_of_phi = 0.0;
    dVdphi = 0.0;
}

#endif /* MYPOTENTIALFUNCTION_HPP */
