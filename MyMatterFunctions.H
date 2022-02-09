/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYMATTERFUNCTIONS_HPP
#define MYMATTERFUNCTIONS_HPP

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

#endif /* MYMATTERFUNCTIONS_HPP */
