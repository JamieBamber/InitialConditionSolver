/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETAIJ_HPP
#define SETAIJ_HPP

#include "PoissonParameters.H"

// Get the bh radius - input loc_bh should be location on the grid relative
// to the centre - note that it gets corrected to the actual location
// (see examples below)
Real get_bh_radius(RealVect &loc_bh, RealVect bh_offset)
{
    loc_bh[0] -= bh_offset[0];
    loc_bh[1] -= bh_offset[1];
    loc_bh[2] -= bh_offset[2];
    return sqrt(loc_bh[0] * loc_bh[0] + loc_bh[1] * loc_bh[1] +
                loc_bh[2] * loc_bh[2]);
}

// Set Aij Bowen York data
// see Alcubierre pg 110 eqn (3.4.22)
Real get_Aij(const int i, const int j, const Real &rbh1, const Real &rbh2,
             const RealVect &n1, const RealVect &n2, const RealVect &J1,
             const RealVect &J2, const RealVect &P1, const RealVect &P2,
             const PoissonParameters &a_params)
{
    std::array<std::array<std::array<double, 3>, 3>, 3> epsilon = {0.};
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[0][2][1] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[1][0][2] = -1.0;

    Real Aij = 1.5 / rbh1 / rbh1 * (n1[i] * P1[j] + n1[j] * P1[i]) +
               1.5 / rbh2 / rbh2 * (n2[i] * P2[j] + n2[j] * P2[i]);
    for (int k = 0; k < SpaceDim; k++)
    {
        Aij +=
            1.5 / rbh1 / rbh1 * (n1[i] * n1[j] - Real(i == j)) * P1[k] * n1[k] +
            1.5 / rbh2 / rbh2 * (n2[i] * n2[j] - Real(i == j)) * P2[k] * n2[k];

        for (int l = 0; l < SpaceDim; l++)
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

void set_binary_bh_Aij(Real Aij[3][3], const IntVect &iv,
                       const RealVect &loc, const PoissonParameters &a_params)
{
    // set the BH values - location
    RealVect loc_bh1 = loc;
    Real rbh1 = get_bh_radius(loc_bh1, a_params.bh1_offset);

    RealVect loc_bh2 = loc;
    Real rbh2 = get_bh_radius(loc_bh2, a_params.bh2_offset);

    RealVect n1 = {loc_bh1[0] / rbh1, loc_bh1[1] / rbh1, loc_bh1[2] / rbh1};
    RealVect n2 = {loc_bh2[0] / rbh2, loc_bh2[1] / rbh2, loc_bh2[2] / rbh2};

    // the Bowen York params
    RealVect J1 = {0.0, 0.0, a_params.bh1_spin};
    RealVect J2 = {0.0, 0.0, a_params.bh2_spin};
    RealVect P1 = {a_params.bh1_momentum[0], a_params.bh1_momentum[1], a_params.bh1_momentum[2]};
    RealVect P2 = {a_params.bh2_momentum[0], a_params.bh2_momentum[1], a_params.bh2_momentum[2]};

    // set the Aij vars
    Aij[0][0] =
        get_Aij(0, 0, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    Aij[1][1] =
        get_Aij(1, 1, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    Aij[2][2] =
        get_Aij(2, 2, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    Aij[0][1] =
        get_Aij(0, 1, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    Aij[0][2] =
        get_Aij(0, 2, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    Aij[1][2] =
        get_Aij(1, 2, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    Aij[1][0] = Aij[0][1];
    Aij[2][0] = Aij[0][2];
    Aij[2][1] = Aij[1][2];
}

// The Brill Lindquist part of psi
Real set_binary_bh_psi(const RealVect &loc, const PoissonParameters &a_params)
{
    // the Bowen York params
    Real m1 = a_params.bh1_bare_mass;
    Real m2 = a_params.bh2_bare_mass;

    // set the BH values - location
    RealVect loc_bh1 = loc;
    Real rbh1 = get_bh_radius(loc_bh1, a_params.bh1_offset);

    RealVect loc_bh2 = loc;
    Real rbh2 = get_bh_radius(loc_bh2, a_params.bh2_offset);

    return m1 / rbh1 + m2 / rbh2;
}

// The part of Aij excluding the Brill Lindquist BH Aij
void set_Aij_reg(Real Aij[3][3], FArrayBox &multigrid_vars_box, const IntVect &iv,
               const RealVect &loc, const RealVect &a_dx,
               const PoissonParameters &a_params,
               const FArrayBox &a_grad_multigrid)
{

    Real trace = (a_grad_multigrid(iv, 3 * c_V0_0 + 0) +
                  a_grad_multigrid(iv, 3 * c_V1_0 + 1) +
                  a_grad_multigrid(iv, 3 * c_V2_0 + 2));

    Aij[0][0] = 0.75 *
            (a_grad_multigrid(iv, 3 * c_V0_0 + 0) +
             a_grad_multigrid(iv, 3 * c_V0_0 + 0)) -
        1. / 2. * trace;

    Aij[1][1] = 0.75 *
            (a_grad_multigrid(iv, 3 * c_V1_0 + 1) +
             a_grad_multigrid(iv, 3 * c_V1_0 + 1)) -
        1. / 2. * trace;

    Aij[2][2] = 0.75 *
            (a_grad_multigrid(iv, 3 * c_V2_0 + 2) +
             a_grad_multigrid(iv, 3 * c_V2_0 + 2)) -
        1. / 2. * trace;

    Aij[0][1] = 0.75 *
                                      (a_grad_multigrid(iv, 3 * c_V0_0 + 1) +
                                       a_grad_multigrid(iv, 3 * c_V1_0 + 0));

    Aij[0][2] = 0.75 *
                                      (a_grad_multigrid(iv, 3 * c_V0_0 + 2) +
                                       a_grad_multigrid(iv, 3 * c_V2_0 + 0));

    Aij[1][2] = 0.75 *
                                      (a_grad_multigrid(iv, 3 * c_V2_0 + 1) +
                                       a_grad_multigrid(iv, 3 * c_V1_0 + 2));
    Aij[1][0] = Aij[0][1];
    Aij[2][0] = Aij[0][2];
    Aij[2][1] = Aij[1][2];
}

/*
void set_deriv_Aij_0(Real d_Aij[3][3][3], FArrayBox &multigrid_vars_box,
                     const IntVect &iv, const RealVect &loc,
                     const RealVect &a_dx, const PoissonParameters &a_params,
                     const FArrayBox &a_grad2_multigrid,
                     const FArrayBox &a_mixed_grad2_multigrid)
{
    Real d2_V[3][3][3];

    d2_V[0][0][0] = a_grad2_multigrid(iv, 3 * c_V0_0 + 0);
    d2_V[0][1][1] = a_grad2_multigrid(iv, 3 * c_V0_0 + 1);
    d2_V[0][2][2] = a_grad2_multigrid(iv, 3 * c_V0_0 + 2);

    d2_V[1][0][0] = a_grad2_multigrid(iv, 3 * c_V1_0 + 0);
    d2_V[1][1][1] = a_grad2_multigrid(iv, 3 * c_V1_0 + 1);
    d2_V[1][2][2] = a_grad2_multigrid(iv, 3 * c_V1_0 + 2);

    d2_V[2][0][0] = a_grad2_multigrid(iv, 3 * c_V2_0 + 0);
    d2_V[2][1][1] = a_grad2_multigrid(iv, 3 * c_V2_0 + 1);
    d2_V[2][2][2] = a_grad2_multigrid(iv, 3 * c_V2_0 + 2);

    d2_V[0][1][0] = a_mixed_grad2_multigrid(iv, 3 * c_V0_0 + 1 + 0 - 1);
    d2_V[0][2][0] = a_mixed_grad2_multigrid(iv, 3 * c_V0_0 + 2 + 0 - 1);
    d2_V[0][1][2] = a_mixed_grad2_multigrid(iv, 3 * c_V0_0 + 1 + 2 - 1);
    d2_V[0][0][1] = d2_V[0][1][0];
    d2_V[0][0][2] = d2_V[0][2][0];
    d2_V[0][2][1] = d2_V[0][1][2];

    d2_V[1][1][0] = a_mixed_grad2_multigrid(iv, 3 * c_V1_0 + 1 + 0 - 1);
    d2_V[1][2][0] = a_mixed_grad2_multigrid(iv, 3 * c_V1_0 + 2 + 0 - 1);
    d2_V[1][1][2] = a_mixed_grad2_multigrid(iv, 3 * c_V1_0 + 1 + 2 - 1);
    d2_V[1][0][1] = d2_V[1][1][0];
    d2_V[1][0][2] = d2_V[1][2][0];
    d2_V[1][2][1] = d2_V[1][1][2];

    d2_V[2][1][0] = a_mixed_grad2_multigrid(iv, 3 * c_V2_0 + 1 + 0 - 1);
    d2_V[2][2][0] = a_mixed_grad2_multigrid(iv, 3 * c_V2_0 + 2 + 0 - 1);
    d2_V[2][1][2] = a_mixed_grad2_multigrid(iv, 3 * c_V2_0 + 1 + 2 - 1);
    d2_V[2][0][1] = d2_V[2][1][0];
    d2_V[2][0][2] = d2_V[2][2][0];
    d2_V[2][2][1] = d2_V[2][1][2];

    Real delta[3][3];
    for (int i = 0; i < SpaceDim; i++)
    {
        for (int j = 0; j < SpaceDim; j++)
        {
            delta[i][j] = 0.0;
        }
    }

    delta[0][0] = 1;
    delta[1][1] = 1;
    delta[2][2] = 1;

    for (int i = 0; i < SpaceDim; i++)
    {
        for (int j = 0; j < SpaceDim; j++)
        {
            for (int k = 0; k < SpaceDim; k++)
            {
                d_Aij[i][j][k] = 0.75 * (d2_V[j][i][k] + d2_V[i][j][k]);

                for (int l = 0; l < SpaceDim; l++)
                {
                    d_Aij[i][j][k] += -0.5 * delta[i][j] * d2_V[l][l][k];
                }
            }
        }
    }
}
*/

#endif /* SETAIJ_HPP */
