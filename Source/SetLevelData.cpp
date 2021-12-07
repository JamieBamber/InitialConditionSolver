/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SetLevelData.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "MyPhiFunction.H"
#include "MyPotentialFunction.H"
#include "PoissonParameters.H"
#include "SetAij.H"
#include "SetLevelDataCpp.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij, scalar field phi,
// its conjugate momentum Pi and the V_i comps which compose Aij
// For now the default setup is 2 Bowen York BHs plus a scalar field
// with some initial user specified configuration
void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars,
                            LevelData<FArrayBox> &a_dpsi, const RealVect &a_dx,
                            const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    const DisjointBoxLayout &grids = a_multigrid_vars.disjointBoxLayout();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &dpsi_box = a_dpsi[dit()];
        for (int comp = 0; comp < NUM_MULTIGRID_VARS; comp++)
        {
            multigrid_vars_box.setVal(0.0, comp);
        }
        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            dpsi_box.setVal(0.0, comp);
        }
        // Iterate over the box
        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            // work out location on the grid
            IntVect iv = bit();
            RealVect loc(iv + 0.5 * RealVect::Unit);
            loc *= a_dx;
            loc -= a_params.domainLength / 2.0;

            // note that we don't include the singular part of psi
            // for the BHs - this is added at the output data stage
            // and when we calculate psi_reg in the rhs etc
            // as it already satisfies Laplacian(psi) = 0
            multigrid_vars_box(iv, c_psi_reg) = a_params.psi_reg;

            // set phi and pi according to user defined function
            multigrid_vars_box(iv, c_phi_0) = my_phi_function(loc, a_params);
            multigrid_vars_box(iv, c_Pi_0) = my_pi_function(loc, a_params);

            // Note that Aij_0 and K_0 are set below since they require
            // gradients
        }
    }
} // end set_initial_conditions

// Set K_0 and check integrands for periodic cases, update Aij
void set_K_and_integrability(LevelData<FArrayBox> &a_integrand,
                             LevelData<FArrayBox> &a_multigrid_vars,
                             const RealVect &a_dx,
                             const PoissonParameters &a_params)
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_integrand.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &integrand_box = a_integrand[dit()];
        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            integrand_box.setVal(0.0, comp);
        }
        Box unghosted_box = integrand_box.box();

        // TODO: Convert these to act at a point not across box
        // and call in bit loop below

        // calculate the laplacian of psi and V_i across the box
        FArrayBox laplace_multigrid(unghosted_box, NUM_CONSTRAINT_VARS);
        get_laplacian(unghosted_box, multigrid_vars_box,
                      Interval(c_psi_reg, c_V3_0), a_dx, laplace_multigrid,
                      a_params);

        // calculate gradients for constructing rho and Aij
        FArrayBox grad_multigrid(unghosted_box, 3 * NUM_MULTIGRID_VARS);
        get_grad(unghosted_box, multigrid_vars_box,
                 Interval(c_psi_reg, c_phi_0), a_dx, grad_multigrid, a_params);

        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc(iv + 0.5 * RealVect::Unit);
            loc *= a_dx;
            loc -= a_params.domainLength / 2.0;

            // integrand = -1.5*term + 1.5 * \bar A_ij \bar A^ij psi_0^-12 +
            // 24 pi rho_grad psi_0^-4  + 12*laplacian(psi_0)*psi^-5
            // JCAurre: Calculate potential and derivative
            Real V_of_phi, dVdphi;
            my_potential_function(V_of_phi, dVdphi,
                                  multigrid_vars_box(iv, c_phi_0), a_params);

            // Calculate the actual value of psi including BH part
            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bh;

            // Assign values of useful quantities
            // TODO: Use Tensor class not arrays of Reals
            Real Pi_0 = multigrid_vars_box(iv, c_Pi_0);
            Real d_phi[3];
            set_deriv1(d_phi, iv, grad_multigrid, c_phi_0);
            Real Aij_reg[3][3];
            Real Aij_bh[3][3];
            set_Aij_reg(Aij_reg, multigrid_vars_box, iv, loc, a_dx, a_params,
                        grad_multigrid);
            set_binary_bh_Aij(Aij_bh, iv, loc, a_params);

            // Compute rhograd from gradients of phi, factors of psi_0
            // accounted for below
            Real rho_gradient = 0;
            for (int i = 0; i < SpaceDim; i++)
            {
                rho_gradient += 0.5 * d_phi[i] * d_phi[i];
            }

            // Also \bar  A_ij \bar A^ij (factors of psi introduced below)
            Real A2_0 = 0.0;
            Real A2_bh = 0.0;
            for (int i = 0; i < SpaceDim; i++)
            {
                for (int j = 0; j < SpaceDim; j++)
                {
                    A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                            (Aij_reg[i][j] + Aij_bh[i][j]);
                    A2_bh += Aij_bh[i][j] * Aij_bh[i][j];
                }
            }

            Real K_0_sq = 0.0;
            if (a_params.is_periodic)
            {
                K_0_sq =
                    1.5 * 8.0 * M_PI * a_params.G_Newton *
                        (pow(Pi_0, 2.0) + 2.0 * V_of_phi) +
                    1.5 * A2_0 * pow(psi_0, -12.0) +
                    24.0 * M_PI * a_params.G_Newton * rho_gradient *
                        pow(psi_0, -4.0) +
                    12.0 * laplace_multigrid(iv, c_psi_reg) * pow(psi_0, -5.0);
            }
            else // in non periodic case exclude psi terms
            {
                K_0_sq =
                    1.5 * 8.0 * M_PI * a_params.G_Newton *
                        (pow(Pi_0, 2.0) + 2.0 * V_of_phi) +
                    // 1.5 * A2_0 * pow(psi_0, -12.0) +
                    24.0 * M_PI * a_params.G_Newton * rho_gradient *
                        pow(psi_0, -4.0);
                    //12.0 * laplace_multigrid(iv, c_psi_reg) * pow(psi_0, -5.0);
            }

            integrand_box(iv, c_psi) =
                -1.5 * (2.0 / 3.0 * K_0_sq -
                        8.0 * M_PI * a_params.G_Newton *
                            (pow(Pi_0, 2.0) + 2.0 * V_of_phi)) +
                1.5 * A2_0 * pow(psi_0, -12.0) +
                24.0 * M_PI * a_params.G_Newton * rho_gradient *
                    pow(psi_0, -4.0) +
                12.0 * laplace_multigrid(iv, c_psi_reg) * pow(psi_0, -5.0);

            integrand_box(iv, c_V1) =
                -8.0 * M_PI * pow(psi_0, 6.0) * Pi_0 * d_phi[0] -
                laplace_multigrid(iv, c_V1_0);

            integrand_box(iv, c_V2) =
                -8.0 * M_PI * pow(psi_0, 6.0) * Pi_0 * d_phi[1] -
                laplace_multigrid(iv, c_V2_0);

            integrand_box(iv, c_V3) =
                -8.0 * M_PI * pow(psi_0, 6.0) * Pi_0 * d_phi[2] -
                laplace_multigrid(iv, c_V3_0);

            // Set value for K
            multigrid_vars_box(iv, c_K_0) = a_params.sign_of_K * sqrt(K_0_sq);
            // be careful when K=0, maybe discontinuity

            // set values for Aij_0
            multigrid_vars_box(iv, c_A11_0) = Aij_reg[0][0] + Aij_bh[0][0];
            multigrid_vars_box(iv, c_A22_0) = Aij_reg[1][1] + Aij_bh[1][1];
            multigrid_vars_box(iv, c_A33_0) = Aij_reg[2][2] + Aij_bh[2][2];
            multigrid_vars_box(iv, c_A12_0) = Aij_reg[0][1] + Aij_bh[0][1];
            multigrid_vars_box(iv, c_A13_0) = Aij_reg[0][2] + Aij_bh[0][2];
            multigrid_vars_box(iv, c_A23_0) = Aij_reg[1][2] + Aij_bh[1][2];
        }
    }
} // end set_K_and_integrability

// set the rhs source for the poisson eqn
void set_rhs(LevelData<FArrayBox> &a_rhs,
             LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx,
             const PoissonParameters &a_params)
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_rhs.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &rhs_box = a_rhs[dit()];
        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            rhs_box.setVal(0.0, comp);
        }
        Box unghosted_box = rhs_box.box();

        // TODO: Make these local functions in bit() and use Tensor class

        // calculate the laplacian of Psi across the box
        FArrayBox laplace_multigrid(unghosted_box, NUM_CONSTRAINT_VARS);
        get_laplacian(unghosted_box, multigrid_vars_box, Interval(c_psi, c_V3),
                      a_dx, laplace_multigrid, a_params);

        // calculate gradients for constructing rho and Aij
        FArrayBox grad_multigrid(unghosted_box, 3 * NUM_MULTIGRID_VARS);
        get_grad(unghosted_box, multigrid_vars_box, Interval(c_psi_reg, c_K_0),
                 a_dx, grad_multigrid, a_params);

        FArrayBox grad2_multigrid(unghosted_box, 3 * NUM_MULTIGRID_VARS);
        get_grad2(unghosted_box, multigrid_vars_box,
                  Interval(c_psi_reg, c_Pi_0), a_dx, grad2_multigrid, a_params);

        FArrayBox mixed_grad2_multigrid(unghosted_box, 3 * NUM_MULTIGRID_VARS);
        get_mixed_grad2(unghosted_box, multigrid_vars_box,
                        Interval(c_psi_reg, c_Pi_0), a_dx,
                        mixed_grad2_multigrid, a_params);

        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc(iv + 0.5 * RealVect::Unit);
            loc *= a_dx;
            loc -= a_params.domainLength / 2.0;

            // rhs = term/8 psi_0^5 - 2 pi rho_grad psi_0  - laplacian(psi_0)
            // JCAurre: Calculate potential and derivative
            Real V_of_phi, dVdphi;
            my_potential_function(V_of_phi, dVdphi,
                                  multigrid_vars_box(iv, c_phi_0), a_params);

            // Calculate useful quantities
            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bh;
            Real K_0 = multigrid_vars_box(iv, c_K_0);
            Real Pi_0 = multigrid_vars_box(iv, c_Pi_0);
            Real d_K[3], d_phi[3], dd_phi[3][3];
            set_deriv1(d_K, iv, grad_multigrid, c_K_0);
            set_deriv1(d_phi, iv, grad_multigrid, c_phi_0);
            set_deriv2(dd_phi, iv, grad2_multigrid, mixed_grad2_multigrid,
                       c_phi_0);

            Real Aij_reg[3][3], Aij_bh[3][3];
            // note that this is \bar A^{ij}, but since conformal metric
            // if flat we have \bar A^{ij} = \bar A_{ij}
            set_Aij_reg(Aij_reg, multigrid_vars_box, iv, loc, a_dx, a_params,
                        grad_multigrid);
            // this is \bar A_{ij}
            set_binary_bh_Aij(Aij_bh, iv, loc, a_params);

            // Compute rhograd from gradients of phi, factors of psi0
            // accounted for below
            Real rho_gradient = 0;
            for (int i = 0; i < SpaceDim; i++)
            {
                rho_gradient += 0.5 * d_phi[i] * d_phi[i];
            }

            // Also \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            for (int i = 0; i < SpaceDim; i++)
            {
                for (int j = 0; j < SpaceDim; j++)
                {
                    A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                            (Aij_reg[i][j] + Aij_bh[i][j]);
                }
            }

            rhs_box(iv, c_psi) =
                0.125 *
                    (2.0 / 3.0 * K_0 * K_0 -
                     8.0 * M_PI * a_params.G_Newton *
                         (pow(Pi_0, 2.0) + 2.0 * V_of_phi)) *
                    pow(psi_0, 5.0) -
                0.125 * A2_0 * pow(psi_0, -7.0) -
                2.0 * M_PI * a_params.G_Newton * rho_gradient * psi_0 -
                laplace_multigrid(iv, c_psi_reg);

            // use the linear form
            rhs_box(iv, c_V1) =
                -8.0 * M_PI * pow(psi_0, 6.0) * Pi_0 * d_phi[0];// -
                //laplace_multigrid(iv, c_V1_0);
            rhs_box(iv, c_V2) =
                -8.0 * M_PI * pow(psi_0, 6.0) * Pi_0 * d_phi[1];// -
                //laplace_multigrid(iv, c_V2_0);
            rhs_box(iv, c_V3) =
                -8.0 * M_PI * pow(psi_0, 6.0) * Pi_0 * d_phi[2];// -
                //laplace_multigrid(iv, c_V3_0);

            rhs_box(iv, c_V1) += 2. / 3. * pow(psi_0, 6.0) * d_K[0];
            rhs_box(iv, c_V2) += 2. / 3. * pow(psi_0, 6.0) * d_K[1];
            rhs_box(iv, c_V3) += 2. / 3. * pow(psi_0, 6.0) * d_K[2];
        }
    }
} // end set_rhs

// set the regrid condition - abs value of this drives AMR regrid
void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                          LevelData<FArrayBox> &a_multigrid_vars,
                          const RealVect &a_dx,
                          const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_condition.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &condition_box = a_condition[dit()];
        condition_box.setVal(0.0, 0);
        Box unghosted_box = condition_box.box();

        // calculate gradients for constructing rho and Aij
        FArrayBox grad_multigrid(unghosted_box, 3 * NUM_MULTIGRID_VARS);
        get_grad(unghosted_box, multigrid_vars_box,
                 Interval(c_psi_reg, c_phi_0), a_dx, grad_multigrid, a_params);

        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc(iv + 0.5 * RealVect::Unit);
            loc *= a_dx;
            loc -= a_params.domainLength / 2.0;

            // Calculate potential and derivative
            Real V_of_phi, dVdphi;
            my_potential_function(V_of_phi, dVdphi,
                                  multigrid_vars_box(iv, c_phi_0), a_params);

            // Calculate useful quantities
            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bh;
            Real K_0 = multigrid_vars_box(iv, c_K_0);
            Real Pi_0 = multigrid_vars_box(iv, c_Pi_0);

            Real d_phi[3];
            set_deriv1(d_phi, iv, grad_multigrid, c_phi_0);
            Real Aij_reg[3][3];
            Real Aij_bh[3][3];
            // note that this is \bar A^{ij}, but since conformal metric
            // if flat we have \bar A^{ij} = \bar A_{ij}
            set_Aij_reg(Aij_reg, multigrid_vars_box, iv, loc, a_dx, a_params,
                        grad_multigrid);
            // this is \bar A_{ij}
            set_binary_bh_Aij(Aij_bh, iv, loc, a_params);

            // JCAurre: compute rhograd from gradients of phi
            Real rho_gradient = 0;
            for (int i = 0; i < SpaceDim; i++)
            {
                rho_gradient += 0.5 * d_phi[i] * d_phi[i];
            }

            // Also \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            for (int i = 0; i < SpaceDim; i++)
            {
                for (int j = 0; j < SpaceDim; j++)
                {
                    A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                            (Aij_reg[i][j] + Aij_bh[i][j]);
                }
            }

            // the condition is similar to the rhs but we take abs
            // value of the contributions and add in BH criteria
            condition_box(iv, 0) =
                1.5 * abs((2.0 / 3.0 * K_0 * K_0 +
                           8.0 * M_PI * a_params.G_Newton *
                               (pow(Pi_0, 2.0) + 2.0 * V_of_phi))) +
                1.5 * A2_0 * pow(psi_0, -7.0) +
                24.0 * M_PI * a_params.G_Newton * abs(rho_gradient) *
                    pow(psi_0, 1.0) + log(psi_0);
            // TODO: Add abs of D_i K and other Mom source terms here
        }
    }
} // end set_regrid_condition

// Add the correction to psi0 after the solver operates
void set_update_psi0(LevelData<FArrayBox> &a_multigrid_vars,
                     LevelData<FArrayBox> &a_dpsi,
                     const Copier &a_exchange_copier)
{
    // first exchange ghost cells for dpsi so they are filled with the correct
    // values
    a_dpsi.exchange(a_dpsi.interval(), a_exchange_copier);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &dpsi_box = a_dpsi[dit()];

        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();

            // Update constraint variables for the linear step
            multigrid_vars_box(iv, c_psi_reg) += dpsi_box(iv, c_psi);
            // these are just linear, so equal instead of add
            multigrid_vars_box(iv, c_V1_0) = dpsi_box(iv, c_V1);
            multigrid_vars_box(iv, c_V2_0) = dpsi_box(iv, c_V2);
            multigrid_vars_box(iv, c_V3_0) = dpsi_box(iv, c_V3);
        }
    }
}

// The coefficient of the I operator on dpsi
void set_a_coef(LevelData<FArrayBox> &a_aCoef,
                LevelData<FArrayBox> &a_multigrid_vars,
                const PoissonParameters &a_params, const RealVect &a_dx)
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &aCoef_box = a_aCoef[dit()];
        for (int iconstraint = 0; iconstraint < NUM_CONSTRAINT_VARS;
             iconstraint++)
        {
            // this prevents small amounts of noise in the sources
            // activating the zero modes - Garfinkle trick!
            // Seems to work best to set this relative to the tolerance
            aCoef_box.setVal(-100.0 * a_params.tolerance, iconstraint);
        }

        // For the non periodic case
        // add back some non trivial psi for the Aij part
        if (!a_params.is_periodic)
        {
            FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
            Box unghosted_box = aCoef_box.box();
            // calculate gradients for constructing rho and Aij
            FArrayBox grad_multigrid(unghosted_box, 3 * NUM_MULTIGRID_VARS);
            get_grad(unghosted_box, multigrid_vars_box,
                     Interval(c_psi_reg, c_phi_0), a_dx, grad_multigrid,
                     a_params);

            BoxIterator bit(unghosted_box);
            for (bit.begin(); bit.ok(); ++bit)
            {
                IntVect iv = bit();
                RealVect loc(iv + 0.5 * RealVect::Unit);
                loc *= a_dx;
                loc -= a_params.domainLength / 2.0;

                // Calculate useful quantities
                Real V_of_phi, dVdphi;
                my_potential_function(V_of_phi, dVdphi,
                                      multigrid_vars_box(iv, c_phi_0),
                                      a_params);

                // Calculate useful quantities
                Real psi_bh = set_binary_bh_psi(loc, a_params);
                Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bh;
                Real K_0 = multigrid_vars_box(iv, c_K_0);
                Real Pi_0 = multigrid_vars_box(iv, c_Pi_0);
                Real d_phi[3];
                set_deriv1(d_phi, iv, grad_multigrid, c_phi_0);

                Real Aij_reg[3][3], Aij_bh[3][3];
                set_Aij_reg(Aij_reg, multigrid_vars_box, iv, loc, a_dx,
                            a_params, grad_multigrid);
                set_binary_bh_Aij(Aij_bh, iv, loc, a_params);

                // Compute rhograd from gradients of phi, factors of psi0
                // accounted for below
                Real rho_gradient = 0;
                for (int i = 0; i < SpaceDim; i++)
                {
                    rho_gradient += 0.5 * d_phi[i] * d_phi[i];
                }
                // Also \bar  A_ij \bar A^ij
                Real A2_0 = 0.0;
                Real A2_bh = 0.0;
                for (int i = 0; i < SpaceDim; i++)
                {
                    for (int j = 0; j < SpaceDim; j++)
                    {
                        A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                                (Aij_reg[i][j] + Aij_bh[i][j]);
                        A2_bh += Aij_bh[i][j] * Aij_bh[i][j];
                    }
                }

                // checked, found errors, should now be right
                aCoef_box(iv, c_psi) =
                    - 0.875 * A2_0 * pow(psi_0, -8.0)
                    - 5.0 / 12.0 * K_0 * K_0 * pow(psi_0, 4.0) +
                    10.0 * M_PI * a_params.G_Newton *
                        (0.5 * Pi_0 * Pi_0 + V_of_phi) * pow(psi_0, 4.0) +
                    10.0 * M_PI * a_params.G_Newton * rho_gradient;
            }
        }
    }
}
// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
// the rhs source of the Poisson eqn
void set_b_coef(LevelData<FArrayBox> &a_bCoef,
                const PoissonParameters &a_params, const RealVect &a_dx)
{

    CH_assert(a_bCoef.nComp() == NUM_CONSTRAINT_VARS);
    int comp_number = 0;

    for (DataIterator dit = a_bCoef.dataIterator(); dit.ok(); ++dit)
    {
        FArrayBox &bCoef_box = a_bCoef[dit()];
        // JCAurre: Loop to set bCoef=1 for all constraint variables
        for (int iconstraint = 0; iconstraint < NUM_CONSTRAINT_VARS;
             iconstraint++)
        {
            bCoef_box.setVal(1.0, iconstraint);
        }
    }
}

// used to set output data for all ADM Vars for GRChombo restart
void set_output_data(LevelData<FArrayBox> &a_grchombo_vars,
                     LevelData<FArrayBox> &a_multigrid_vars,
                     const PoissonParameters &a_params, const RealVect &a_dx)
{

    CH_assert(a_grchombo_vars.nComp() == NUM_GRCHOMBO_VARS);
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    const DisjointBoxLayout &grids = a_grchombo_vars.disjointBoxLayout();
    DataIterator dit = a_grchombo_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &grchombo_vars_box = a_grchombo_vars[dit()];
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // first set everything to zero
        for (int comp = 0; comp < NUM_GRCHOMBO_VARS; comp++)
        {
            grchombo_vars_box.setVal(0.0, comp);
        }

        // now set non zero terms - const across whole box
        // Conformally flat, and lapse = 1
        grchombo_vars_box.setVal(1.0, c_h11);
        grchombo_vars_box.setVal(1.0, c_h22);
        grchombo_vars_box.setVal(1.0, c_h33);
        grchombo_vars_box.setVal(1.0, c_lapse);

        // now non constant terms by location
        Box this_box = grchombo_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc(iv + 0.5 * RealVect::Unit);
            loc *= a_dx;
            loc -= a_params.domainLength / 2.0;

            // GRChombo conformal factor chi = psi^-4
            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real chi = pow(multigrid_vars_box(iv, c_psi_reg) + psi_bh, -4.0);
            grchombo_vars_box(iv, c_chi) = chi;
            Real factor = pow(chi, 1.5);

            // Copy phi and Aij across - note this is now \tilde Aij not
            // \bar Aij
            grchombo_vars_box(iv, c_phi) = multigrid_vars_box(iv, c_phi_0);
            grchombo_vars_box(iv, c_Pi) = multigrid_vars_box(iv, c_Pi_0);

            grchombo_vars_box(iv, c_K) = multigrid_vars_box(iv, c_K_0);
            grchombo_vars_box(iv, c_A11) =
                multigrid_vars_box(iv, c_A11_0) * factor;
            grchombo_vars_box(iv, c_A12) =
                multigrid_vars_box(iv, c_A12_0) * factor;
            grchombo_vars_box(iv, c_A13) =
                multigrid_vars_box(iv, c_A13_0) * factor;
            grchombo_vars_box(iv, c_A22) =
                multigrid_vars_box(iv, c_A22_0) * factor;
            grchombo_vars_box(iv, c_A23) =
                multigrid_vars_box(iv, c_A23_0) * factor;
            grchombo_vars_box(iv, c_A33) =
                multigrid_vars_box(iv, c_A33_0) * factor;
        }
    }
}
