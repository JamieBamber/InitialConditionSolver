/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYMATTERFUNCTIONS_HPP
#define MYMATTERFUNCTIONS_HPP

Real my_phi_function(const RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  if(!a_params.is_periodic)
  {
    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return a_params.phi_0 +
         a_params.phi_amplitude *
             exp(-0.5 * pow(rr / a_params.phi_wavelength, 2.0));
  }
  else
  // sine waves in 3 directions - use with periodic BCs
  {
      return a_params.phi_0 +
          a_params.phi_amplitude / 3.0 *
              (cos(2 * M_PI * loc[0] * a_params.phi_wavelength / L[0]) +
               cos(2 * M_PI * loc[1] * a_params.phi_wavelength / L[1]) +
               cos(2 * M_PI * loc[2] * a_params.phi_wavelength / L[2]));
  }
}

Real my_Pi_function(const RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  if(!a_params.is_periodic)
  {
    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return a_params.pi_0 + a_params.pi_amplitude *
                             exp(-0.5 * pow(rr / a_params.pi_wavelength, 2.0));
  }
  else
  {
    // sine waves in 3 directions - use with periodic BCs
    return a_params.pi_0 +
          a_params.pi_amplitude / 3.0 *
              (cos(2 * M_PI * loc[0] * a_params.pi_wavelength / L[0]) +
               cos(2 * M_PI * loc[1] * a_params.pi_wavelength / L[1]) +
               cos(2 * M_PI * loc[2] * a_params.pi_wavelength / L[2]));
  }
}

Real my_potential_function(const Real &phi_here,
                           const PoissonParameters &a_params)
{
    Real mpl = 1.0;
    return a_params.pot_Lambda * 
               pow(1. - exp((phi_here / a_params.pot_mu * pow(mpl, -1.0))), 2.0);
}

#endif /* MYMATTERFUNCTIONS_HPP */
