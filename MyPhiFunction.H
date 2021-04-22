/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYPHIFUNCTION_HPP
#define MYPHIFUNCTION_HPP

Real my_phi_function(RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
  return a_params.phi_0 +
         a_params.phi_amplitude *
             exp(-0.5 * pow(rr / a_params.phi_wavelength, 2.0));

  // sine waves in 3 directions - use with periodic BCs
  // return a_params.phi_0 +
  //        a_params.phi_amplitude / 3.0 *
  //            (cos(2 * M_PI * loc[0] * a_params.phi_wavelength / L[0]) +
  //             cos(2 * M_PI * loc[1] * a_params.phi_wavelength / L[1]) +
  //             cos(2 * M_PI * loc[2] * a_params.phi_wavelength / L[2]));
}

Real my_pi_function(RealVect loc, const PoissonParameters &a_params) {

  RealVect L = a_params.domainLength;

  // gaussian - use with Dirichlet BCs
  Real rr = sqrt(D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]));
  return a_params.pi_0 + a_params.pi_amplitude *
                             exp(-0.5 * pow(rr / a_params.pi_wavelength, 2.0));

  // sine waves in 3 directions - use with periodic BCs
  // return a_params.pi_0 +
  //        a_params.pi_amplitude / 3.0 *
  //            (cos(2 * M_PI * loc[0] * a_params.pi_wavelength / L[0]) +
  //             cos(2 * M_PI * loc[1] * a_params.pi_wavelength / L[1]) +
  //             cos(2 * M_PI * loc[2] * a_params.pi_wavelength / L[2]));
}

#endif /* MYPHIFUNCTION_HPP */
