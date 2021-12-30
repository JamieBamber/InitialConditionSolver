/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETPSIANDAIJ_HPP
#define SETPSIANDAIJ_HPP

#include "PoissonParameters.H"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "SetLevelDataGradients.H"

// Get the bh radius - input loc_bh should be location on the grid relative
// to the centre - note that it gets corrected to the actual location
// (see examples below)
void get_bh_coords(Real &bh_radius, RealVect &loc_bh, 
                   const RealVect &loc, const Real bh_x_offset)
{
    // set coords
    FOR1(i) {loc_bh[i] = loc[i];}
    loc_bh[0] -= bh_x_offset;

    // set radius
    Real bh_radius_squared = 0.0;
    FOR1(i) {bh_radius_squared += loc_bh[i] * loc_bh[i];}
    bh_radius = sqrt(bh_radius_squared);
}

// The Brill Lindquist part of psi
Real set_binary_bh_psi(const RealVect &loc, const PoissonParameters &a_params)
{
    // the Bowen York params
    Real m1 = a_params.bh1_bare_mass;
    Real m2 = a_params.bh2_bare_mass;

    // set the BH values - location
    RealVect loc_bh1;
    Real rbh1;
    get_bh_coords(rbh1, loc_bh1, loc, a_params.bh1_offset);

    RealVect loc_bh2;
    Real rbh2;
    get_bh_coords(rbh2, loc_bh2, loc, a_params.bh2_offset);

    return m1 / rbh1 + m2 / rbh2;
}

// Set Aij Bowen York data
// see Alcubierre pg 110 eqn (3.4.22)
Real get_Aij_component(const int i, const int j,
             const Real &rbh1, const Real &rbh2,
             const RealVect &n1, const RealVect &n2, const RealVect &J1,
             const RealVect &J2, const RealVect &P1, const RealVect &P2,
             const PoissonParameters &a_params)
{
    Tensor<3, Real> epsilon = TensorAlgebra::epsilon();

    Real Aij = 1.5 / rbh1 / rbh1 * (n1[i] * P1[j] + n1[j] * P1[i]) +
               1.5 / rbh2 / rbh2 * (n2[i] * P2[j] + n2[j] * P2[i]);
    FOR1(k)
    {
        Aij +=
            1.5 / rbh1 / rbh1 * (n1[i] * n1[j] - Real(i == j)) * P1[k] * n1[k] +
            1.5 / rbh2 / rbh2 * (n2[i] * n2[j] - Real(i == j)) * P2[k] * n2[k];

        FOR1(l)
        {
            Aij += -3.0 / rbh1 / rbh1 / rbh1 *
                       (epsilon[i][l][k] * n1[j] + epsilon[j][l][k] * n1[i]) *
                       n1[l] * J1[k] -
                   3.0 / rbh2 / rbh2 / rbh2 *
                       (epsilon[i][l][k] * n2[j] + epsilon[j][l][k] * n2[i]) *
                       n2[l] * J2[k];
        }
    }
    return Aij;
}

void set_Aij_bh(Tensor<2, Real> Aij,
                const RealVect &loc, const PoissonParameters &a_params)
{
    // get the BH values - location and radius
    RealVect loc_bh1;
    Real rbh1;
    get_bh_coords(rbh1, loc_bh1, loc, a_params.bh1_offset);

    RealVect loc_bh2;
    Real rbh2;
    get_bh_coords(rbh2, loc_bh2, loc, a_params.bh2_offset);

    // normal vectors
    RealVect n1, n2;
    FOR1(i)
    {
        n1[i] = loc_bh1[i] / rbh1;
        n2[i] = loc_bh2[i] / rbh2;
    }

    // the Bowen York params
    // KC TODO: Make these general vectors
    RealVect J1 = {0.0, 0.0, a_params.bh1_spin};
    RealVect J2 = {0.0, 0.0, a_params.bh2_spin};
    RealVect P1 = {0.0, a_params.bh1_momentum, 0.0};
    RealVect P2 = {0.0, a_params.bh2_momentum, 0.0};

    // set the Aij vars
    FOR2(i,j)
    {
        Aij[i][j] =
            get_Aij_component(i, j, rbh1, rbh2, n1, n2, 
                                  J1, J2, P1, P2, a_params);
    }
}

// The part of Aij excluding the Brill Lindquist BH Aij
void set_Aij_reg(Tensor<2, Real> Aij, 
                 FArrayBox &multigrid_vars_box, 
                 const IntVect &iv,
                 const RealVect &a_dx,
                 const PoissonParameters &a_params)
{
    // get the derivs
    Tensor<2, Real, SpaceDim> d2_U = 
        get_d2(iv, multigrid_vars_box, a_dx, c_U_0);
    Tensor<1, Real, SpaceDim> d1_V1 = 
        get_d1(iv, multigrid_vars_box, a_dx, c_V1_0);
    Tensor<1, Real, SpaceDim> d1_V2 = 
        get_d1(iv, multigrid_vars_box, a_dx, c_V2_0);
    Tensor<1, Real, SpaceDim> d1_V3 = 
        get_d1(iv, multigrid_vars_box, a_dx, c_V3_0);

    using namespace TensorAlgebra;
    Real trace = 0.0;
    FOR1(i)
    {
        trace = delta(i,0) * d1_V1[i] 
              + delta(i,1) * d1_V2[i]
              + delta(i,2) * d1_V3[i];
              + d2_U[i][i];
    }

    // note derivative index j is the last one
    Tensor<2, Real, SpaceDim> d1_Vi;
    FOR2(i,j)
    {
        d1_Vi[i][j] = delta(i,0) * d1_V1[j]
                    + delta(i,1) * d1_V2[j]
                    + delta(i,2) * d1_V3[j];
    }

    // set the values of Aij
    FOR2(i,j)
    {
        Aij[i][j] = d1_Vi[i][j] + d2_U[i][j]
                  + d1_Vi[j][i] + d2_U[j][i]
                  - 2.0 / 3.0 * delta(i,j) * trace;
    }
}

#endif /* SETPSIANDAIJ_HPP */
