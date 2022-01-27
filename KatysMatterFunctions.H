/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYMATTERFUNCTIONS_HPP
#define MYMATTERFUNCTIONS_HPP

Real my_phi_Re_function(const RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  if(!a_params.periodic_directions_exist)
  {
    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return a_params.phi_0 +
         a_params.phi_amplitude * loc[0] * rr *
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

Real my_phi_Im_function(const RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  if(!a_params.periodic_directions_exist)
  {
    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return 0;
  }
  else
  // sine waves in 3 directions - use with periodic BCs
  {
      return a_params.phi_0 +
          a_params.phi_amplitude / 3.0 *
              (sin(2 * M_PI * loc[0] * a_params.phi_wavelength / L[0]) +
               sin(2 * M_PI * loc[1] * a_params.phi_wavelength / L[1]) +
               sin(2 * M_PI * loc[2] * a_params.phi_wavelength / L[2]));
  }
}

Real my_Pi_Re_function(const RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  if(!a_params.periodic_directions_exist)
  {
    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return a_params.pi_0 - a_params.pi_amplitude * loc[1] * rr *
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

Real my_Pi_Im_function(const RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  if(!a_params.periodic_directions_exist)
  {
    Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
    return 0.0;
  }
  else
  {
    // sine waves in 3 directions - use with periodic BCs
    return a_params.pi_0 +
          a_params.pi_amplitude / 3.0 *
              (sin(2 * M_PI * loc[0] * a_params.pi_wavelength / L[0]) +
               sin(2 * M_PI * loc[1] * a_params.pi_wavelength / L[1]) +
               sin(2 * M_PI * loc[2] * a_params.pi_wavelength / L[2]));
  }
}

Real my_potential_function(const Real &phi_Re_here, const Real &phi_Im_here,
                           const PoissonParameters &a_params)
{
    return a_params.pot_mu * a_params.pot_mu * (phi_Re_here * phi_Re_here + phi_Im_here * phi_Im_here);
}

#endif /* MYMATTERFUNCTIONS_HPP */
