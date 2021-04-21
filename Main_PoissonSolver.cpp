/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "DebugDump.H"
#include "FABView.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "MultilevelLinearOp.H"
#include "ParmParse.H"
#include "PoissonParameters.H"
#include "SetBCs.H"
#include "SetGrids.H"
#include "SetLevelData.H"
#include "UsingNamespace.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "WriteOutput.H"
#include "computeNorm.H"
#include "computeSum.H"
#include <iostream>

#ifdef CH_Linux
// Should be undefined by default
//#define TRAP_FPE
#undef TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

using std::cerr;

// Sets up and runs the solver
// The equation solved is: [aCoef*I + bCoef*Laplacian](dpsi) = rhs
// We assume conformal flatness, K=const and Momentum constraint satisfied
// by chosen Aij (for now sourced by Bowen York data),
// lapse = 1 shift = 0, phi is the scalar field and is used to
// calculate the rhs, Pi = dphidt = 0.
int poissonSolve(const Vector<DisjointBoxLayout> &a_grids,
                 const PoissonParameters &a_params)
{

    // the params reader
    ParmParse pp;

    // create the necessary hierarchy of data components
    int nlevels = a_params.numLevels;
    // the user set initial conditions - currently including psi, phi, A_ij
    Vector<LevelData<FArrayBox> *> multigrid_vars(nlevels, NULL);
    // the correction to the conformal factor - what the solver solves for
    Vector<LevelData<FArrayBox> *> dpsi(nlevels, NULL);
    // the solver vars - coefficients and source
    Vector<LevelData<FArrayBox> *> rhs(nlevels, NULL);
    // the integrand for constant K integrability condition
    Vector<LevelData<FArrayBox> *> integrand(nlevels, NULL);
    // the coeff for the I term
    Vector<RefCountedPtr<LevelData<FArrayBox>>> aCoef(nlevels);
    // the coeff for the Laplacian
    Vector<RefCountedPtr<LevelData<FArrayBox>>> bCoef(nlevels);

    // Grid params
    Vector<ProblemDomain> vectDomains(nlevels); // the domains
    Vector<RealVect> vectDx(nlevels); // the grid spacings on each level

    // Set temp vars at coarsest level values
    RealVect dxLev = RealVect::Unit;
    dxLev *= a_params.coarsestDx;
    ProblemDomain domLev(a_params.coarsestDomain);
    IntVect ghosts = 3 * IntVect::Unit;

    // Declare variables here, with num comps, and ghosts for all
    // sources NB - we want output data to have 3 ghost cells to match GRChombo,
    // although not currently needed for 2nd order stencils used here
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
        multigrid_vars[ilev] =
            new LevelData<FArrayBox>(a_grids[ilev], NUM_MULTIGRID_VARS, ghosts);
        dpsi[ilev] = new LevelData<FArrayBox>(a_grids[ilev],
                                              NUM_CONSTRAINTS_VARS, ghosts);
        rhs[ilev] = new LevelData<FArrayBox>(
            a_grids[ilev], NUM_CONSTRAINTS_VARS, IntVect::Zero);
        integrand[ilev] = new LevelData<FArrayBox>(
            a_grids[ilev], NUM_CONSTRAINTS_VARS, IntVect::Zero);
        aCoef[ilev] =
            RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(
                a_grids[ilev], NUM_CONSTRAINTS_VARS, IntVect::Zero));
        bCoef[ilev] =
            RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(
                a_grids[ilev], NUM_CONSTRAINTS_VARS, IntVect::Zero));
        vectDomains[ilev] = domLev;
        vectDx[ilev] = dxLev;
        // set initial guess for psi and zero dpsi
        // and values for other multigrid sources - phi and Aij
        set_initial_conditions(*multigrid_vars[ilev], *dpsi[ilev], vectDx[ilev],
                               a_params);

        // prepare temp dx, domain vars for next level
        dxLev /= a_params.refRatio[ilev];
        domLev.refine(a_params.refRatio[ilev]);
    }

    // set up linear operator
    int lBase = 0;
    MultilevelLinearOp<FArrayBox> mlOp;
    BiCGStabSolver<Vector<LevelData<FArrayBox> *>>
        solver; // define solver object

    // default or read in solver params
    int numMGIter = 1;
    pp.query("numMGIterations", numMGIter);
    mlOp.m_num_mg_iterations = numMGIter;

    int numMGSmooth = 4;
    pp.query("numMGsmooth", numMGSmooth);
    mlOp.m_num_mg_smooth = numMGSmooth;

    int preCondSolverDepth = -1;
    pp.query("preCondSolverDepth", preCondSolverDepth);
    mlOp.m_preCondSolverDepth = preCondSolverDepth;

    Real tolerance = 1.0e-7;
    pp.query("tolerance", tolerance);

    int max_iter = 10;
    pp.query("max_iterations", max_iter);

    int max_NL_iter = 4;
    pp.query("max_NL_iterations", max_NL_iter);

    // Iterate linearised Poisson eqn for NL solution
    Real dpsi_norm0 = 0.0;
	Real dpsi_norm1 = 1.0;
    Real constant_K = 0.0;
    for (int NL_iter = 0; NL_iter < max_NL_iter; NL_iter++)
    {

        pout() << "Main Loop Iteration " << (NL_iter + 1) << " out of "
               << max_NL_iter << endl;

        // Set integrability condition on K if periodic
        if (a_params.periodic[0] == 1)
        {
            // Calculate values for integrand here with K unset
            pout() << "Computing average K value... " << endl;
            for (int ilev = 0; ilev < nlevels; ilev++)
            {
                set_constant_K_integrand(*integrand[ilev],
                                         *multigrid_vars[ilev], vectDx[ilev],
                                         a_params);
            }
            Real integral_S1 = computeSum(integrand, a_params.refRatio,
                                       a_params.coarsestDx, Interval(1, 1));
            Real integral_S2 = computeSum(integrand, a_params.refRatio,
                                       a_params.coarsestDx, Interval(2, 2));
            Real integral_S3 = computeSum(integrand, a_params.refRatio,
                                       a_params.coarsestDx, Interval(3, 3));
            Real volume = a_params.domainLength[0] * a_params.domainLength[1] *
                          a_params.domainLength[2];
            constant_K = 0.0;//M_PI * abs(integral) / volume;
            pout() << "Integral of S1 " << integral_S1 << endl;
            pout() << "Integral of S2 " << integral_S2 << endl;
            pout() << "Integral of S3 " << integral_S3 << endl;
        }

        // Calculate values for coefficients here - see SetLevelData.cpp
        // for details
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            set_a_coef(*aCoef[ilev], *multigrid_vars[ilev], a_params,
                       vectDx[ilev], constant_K);
            set_b_coef(*bCoef[ilev], a_params, vectDx[ilev]);
            set_rhs(*rhs[ilev], *multigrid_vars[ilev], vectDx[ilev], a_params,
                    constant_K);
        }
        // set up solver factory
        RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox>>> opFactory =
            RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox>>>(
                defineOperatorFactory(a_grids, vectDomains, aCoef, bCoef,
                                      a_params));

        // define the multi level operator
        mlOp.define(a_grids, a_params.refRatio, vectDomains, vectDx, opFactory,
                    lBase);
        // set the more solver params
        bool homogeneousBC = false;
        solver.define(&mlOp, homogeneousBC);
        solver.m_verbosity = a_params.verbosity;
        solver.m_normType = 0;
        solver.m_eps = tolerance;
        solver.m_imax = max_iter;

        // output the data before the solver acts to check starting conditions
        output_solver_data(dpsi, rhs, multigrid_vars, a_grids, a_params,
                           NL_iter);

        // Engage!
        solver.solve(dpsi, rhs);

        // Add the solution to the linearised eqn to the previous iteration
        // ie psi -> psi + dpsi
        // need to fill interlevel and intralevel ghosts first in dpsi
        for (int ilev = 0; ilev < nlevels; ilev++)
        {

            // For interlevel ghosts
            if (ilev > 0)
            {
                QuadCFInterp quadCFI(a_grids[ilev], &a_grids[ilev - 1],
                                     vectDx[ilev][0], a_params.refRatio[ilev],
                                     NUM_CONSTRAINTS_VARS, vectDomains[ilev]);
                quadCFI.coarseFineInterp(*dpsi[ilev], *dpsi[ilev - 1]);
            }

		}
		Real max_dpsi = computeMax(dpsi, a_params.refRatio, Interval(0, 0));
		pout() << "max dpsi" << max_dpsi << endl;
		for (int ilev = 0; ilev < nlevels; ilev++)
		{
            // For intralevel ghosts - this is done in set_update_phi0
            // but need the exchange copier object to do this
            Copier exchange_copier;
            exchange_copier.exchangeDefine(a_grids[ilev], ghosts);

            // now the update
            set_update_psi0(*multigrid_vars[ilev], *dpsi[ilev],
                            exchange_copier, 0*max_dpsi);

// removed for now
//            set_constant_K_integrand(*integrand[ilev],
//                                         *multigrid_vars[ilev], vectDx[ilev],
//                                         a_params);
//            set_rhs(*rhs[ilev], *multigrid_vars[ilev], vectDx[ilev], a_params,
//                    constant_K);
        }

        // check if converged or diverging and if so exit NL iteration for loop
        dpsi_norm0 =
            computeNorm(rhs, a_params.refRatio, a_params.coarsestDx,
                        Interval(0, 0)); // TODO JCAurre: not completely sure
        dpsi_norm1 =
            computeNorm(rhs, a_params.refRatio, a_params.coarsestDx,
                        Interval(1, 3)); // TODO JCAurre: not completely sure
        pout() << "The norm of rhs Ham after step " << NL_iter + 1 << " is "
               << dpsi_norm0 << endl;
        pout() << "The norm of rhs Mom after step " << NL_iter + 1 << " is "
               << dpsi_norm1 << endl;
 
		if ((dpsi_norm0 < tolerance  && dpsi_norm1 < tolerance ) || dpsi_norm0 > 1e5 || dpsi_norm1 > 1e5)
        {
            break;
        }

    } // end NL iteration loop

    pout() << "The norm of rhs Ham at the final step was " << dpsi_norm0 << endl;
    pout() << "The norm of rhs Mom at the final step was " << dpsi_norm1 << endl;

	// Mayday if result not converged at all - using a fairly generous threshold
    // for this as usually non convergence means everything goes nuts
    if (dpsi_norm0 > 1e-1 || dpsi_norm1 > 1e-1)
    {
        MayDay::Error(
            "NL iterations did not converge - may need a better initial guess");
    }

    // now output final data in a form which can be read as a checkpoint file
    // for the GRChombo AMR time dependent runs
    output_final_data(multigrid_vars, a_grids, vectDx, vectDomains, a_params,
                      constant_K);

    // clean up data
    for (int level = 0; level < multigrid_vars.size(); level++)
    {
        if (multigrid_vars[level] != NULL)
        {
            delete multigrid_vars[level];
            multigrid_vars[level] = NULL;
        }
        if (rhs[level] != NULL)
        {
            delete rhs[level];
            rhs[level] = NULL;
        }
        if (integrand[level] != NULL)
        {
            delete integrand[level];
            integrand[level] = NULL;
        }
        if (dpsi[level] != NULL)
        {
            delete dpsi[level];
            dpsi[level] = NULL;
        }
    }

    int exitStatus = solver.m_exitStatus;
    // note that for AMRMultiGrid, success = 1.
    exitStatus -= 1;
    return exitStatus;
}

// Main function - keep this simple with just setup and read params
int main(int argc, char *argv[])
{
    int status = 0;
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif
    // Scoping trick
    {
        if (argc < 2)
        {
            cerr << " usage " << argv[0] << " <input_file_name> " << endl;
            exit(0);
        }

        char *inFile = argv[1];
        ParmParse pp(argc - 2, argv + 2, NULL, inFile);

        PoissonParameters params;
        Vector<DisjointBoxLayout> grids;

        // read params from file
        getPoissonParameters(params);

        // set up the grids, using the rhs for tagging to decide
        // where needs additional levels
        set_grids(grids, params);

        // Solve the equations!
        status = poissonSolve(grids, params);

    } // End scoping trick

#ifdef CH_MPI
    MPI_Finalize();
#endif
    return status;
}
