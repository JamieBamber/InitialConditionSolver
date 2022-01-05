#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SetBCs.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Global BCRS definitions
std::vector<int> GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int> GlobalBCRS::s_bcHi = std::vector<int>();
bool GlobalBCRS::s_areBCsParsed = false;
bool GlobalBCRS::s_valueParsed_single = false;
bool GlobalBCRS::s_valueParsed_vector = false;

// BCValueHolder functions, ie a void-type function with the 4
// arguments given pos [x,y,z] position on center of cell edge int dir
// direction, x being 0 int side -1 for low, +1 = high, fill in the a_values
// array
// For the moment these assume the value is constant for all positions and dirs
void ParseSingleValue(Real *pos, int *dir, Side::LoHiSide *side, Real *a_values)
{
    ParmParse pp;
    Real bcVal = 0.0;
    pp.query("bc_value", bcVal);
    a_values[0] = bcVal;
    GlobalBCRS::s_valueParsed_single = true;
}

void ParseVectorValue(Real *pos, int *dir, Side::LoHiSide *side, Real *a_values)
{
    ParmParse pp;
    Real bcVal_vector = 0.0;
    pp.query("bc_value_vector", bcVal_vector);
    for (int icomp = 0; icomp < NUM_CONSTRAINT_VARS; icomp++)
    {
        a_values[icomp] = bcVal_vector;
    }
    GlobalBCRS::s_valueParsed_vector = true;
}

void ParseBC(FArrayBox &a_state, const Box &a_valid,
             const ProblemDomain &a_domain, Real a_dx, bool a_homogeneous)
{
    if (!a_domain.domainBox().contains(a_state.box()))
    {
        if (!GlobalBCRS::s_areBCsParsed)
        {
            ParmParse pp;
            pp.getarr("lo_boundary", GlobalBCRS::s_bcLo, 0, SpaceDim);
            pp.getarr("hi_boundary", GlobalBCRS::s_bcHi, 0, SpaceDim);
            GlobalBCRS::s_areBCsParsed = true;
        }

        const Box valid = a_valid;
        BCValueHolder single_bc(ParseSingleValue); // pointer to void function
        BCValueHolder vector_bc(ParseVectorValue); // pointer to void function

        // this boundary condition leaves the boundary cells set as they are
        // initially (so zeros as set in the InitialConditions)
        // doNothingBC(a_state, valid, a_domain.domainBox(), a_dx,
        // a_homogeneous);

        // This sets all the comps to zero in the BCs
        // which is the default
        ZeroBC(a_state, valid, a_domain.domainBox(), a_dx, a_homogeneous,
               a_state.interval());

        // One can impose alternatives
        // on different components, using the functions defined in BCFunc.H
        // Here we will use them to impose reflective BCs by setting even parity
        // vars to Neumann and odd parity vars to Dirichlet
        for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
            // not periodic and reflective
            if (!a_domain.isPeriodic(idir) && (GlobalBCRS::s_bcLo[idir] == 1 ||
                                               GlobalBCRS::s_bcHi[idir] == 1))
            {
                int num_ghosts = 1;
                Box ghostBoxLo = adjCellBox(valid, idir, Side::Lo, num_ghosts);
                Box ghostBoxHi = adjCellBox(valid, idir, Side::Hi, num_ghosts);
                if (!a_domain.domainBox().contains(ghostBoxLo) &&
                    GlobalBCRS::s_bcLo[idir] == 1)
                {
                    Interval psi_comp = Interval(c_psi, c_psi);
                    NeumBC(a_state, valid, a_dx, a_homogeneous, single_bc, idir,
                           Side::Lo, psi_comp);

                    FOR1(Vi_dir)
                    {
                        Interval Vi_comp =
                            Interval(c_V1 + Vi_dir, c_V1 + Vi_dir);
                        if (Vi_dir == idir)
                        {
                            DiriBC(a_state, valid, a_dx, a_homogeneous,
                                   single_bc, idir, Side::Lo, Vi_comp);
                        }
                        else
                        {
                            NeumBC(a_state, valid, a_dx, a_homogeneous,
                                   single_bc, idir, Side::Lo, Vi_comp);
                        }
                    }

                    Interval U_comp = Interval(c_U, c_U);
                    NeumBC(a_state, valid, a_dx, a_homogeneous, single_bc, idir,
                           Side::Lo, U_comp);
                }

                if (!a_domain.domainBox().contains(ghostBoxHi) &&
                    GlobalBCRS::s_bcLo[idir] == 1)
                {
                    Interval psi_comp = Interval(c_psi, c_psi);
                    NeumBC(a_state, valid, a_dx, a_homogeneous, single_bc, idir,
                           Side::Hi, psi_comp);

                    FOR1(Vi_dir)
                    {
                        Interval Vi_comp =
                            Interval(c_V1 + Vi_dir, c_V1 + Vi_dir);
                        if (Vi_dir == idir)
                        {
                            DiriBC(a_state, valid, a_dx, a_homogeneous,
                                   single_bc, idir, Side::Hi, Vi_comp);
                        }
                        else
                        {
                            NeumBC(a_state, valid, a_dx, a_homogeneous,
                                   single_bc, idir, Side::Hi, Vi_comp);
                        }
                    }

                    Interval U_comp = Interval(c_U, c_U);
                    NeumBC(a_state, valid, a_dx, a_homogeneous, single_bc, idir,
                           Side::Hi, U_comp);
                }
            } // else - is periodic
        }     // close for idir
    }
}

void ZeroBC(FArrayBox &a_state, const Box &a_valid,
            const ProblemDomain &a_domain, Real a_dx, bool a_homogeneous,
            Interval a_interval)
{
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        for (SideIterator sit; sit.ok(); ++sit)
        {
            int num_ghosts = 1;
            Box toRegion = adjCellBox(a_valid, idir, sit(), num_ghosts);
            toRegion &= a_state.box();

            for (BoxIterator bit(toRegion); bit.ok(); ++bit)
            {
                const IntVect &ivTo = bit();

                for (int icomp = a_interval.begin(); icomp <= a_interval.end();
                     icomp++)
                {
                    // Set BCs to zero
                    a_state(ivTo, icomp) = 0.0;
                }
            }
        }
    }
}
