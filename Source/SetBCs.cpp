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
std::vector<bool> GlobalBCRS::s_printedThatLo_psi =
    std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi_psi =
    std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatLo_Vi =
    std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi_Vi =
    std::vector<bool>(SpaceDim, false);
std::vector<int> GlobalBCRS::s_bcLo_psi = std::vector<int>();
std::vector<int> GlobalBCRS::s_bcHi_psi = std::vector<int>();
std::vector<int> GlobalBCRS::s_bcLo_Vi = std::vector<int>();
std::vector<int> GlobalBCRS::s_bcHi_Vi = std::vector<int>();
RealVect GlobalBCRS::s_trigvec = RealVect::Zero;
bool GlobalBCRS::s_areBCsParsed_psi = false;
bool GlobalBCRS::s_valueParsed_psi = false;
bool GlobalBCRS::s_areBCsParsed_Vi = false;
bool GlobalBCRS::s_valueParsed_Vi = false;
bool GlobalBCRS::s_trigParsed = false;

// BCValueHolder class, which is a pointer to a void-type function with the 4
// arguements given pos [x,y,z] position on center of cell edge int dir
// direction, x being 0 int side -1 for low, +1 = high, fill in the a_values
// array

void ParseValuePsi(Real *pos, int *dir, Side::LoHiSide *side, Real *a_values)
{
    ParmParse pp;
    Real bcVal_psi;
    pp.get("bc_value_psi", bcVal_psi);
    a_values[0] = bcVal_psi;
}

void ParseValueVi(Real *pos, int *dir, Side::LoHiSide *side, Real *a_values)
{
    ParmParse pp;
    Real bcVal_Vi;
    pp.get("bc_value_Vi", bcVal_Vi);
    a_values[0] = bcVal_Vi;
    a_values[1] = bcVal_Vi;
    a_values[2] = bcVal_Vi;
}

void ParseBC(FArrayBox &a_state, const Box &a_valid,
             const ProblemDomain &a_domain, Real a_dx, bool a_homogeneous)
{
    if (!a_domain.domainBox().contains(a_state.box()))
    {

        if (!GlobalBCRS::s_areBCsParsed)
        {
            ParmParse pp;
            pp.getarr("bc_lo_psi", GlobalBCRS::s_bcLo_psi, 0, SpaceDim);
            pp.getarr("bc_hi_psi", GlobalBCRS::s_bcHi_psi, 0, SpaceDim);
            pp.getarr("bc_lo_Vi", GlobalBCRS::s_bcLo_Vi, 0, SpaceDim);
            pp.getarr("bc_hi_Vi", GlobalBCRS::s_bcHi_Vi, 0, SpaceDim);
            GlobalBCRS::s_areBCsParsed = true;
        }

        Box valid = a_valid;

        for (int i = 0; i < CH_SPACEDIM; ++i)
        {

            // periodic? If not, check if Dirichlet or Neumann
            if (!a_domain.isPeriodic(i))
            {
                Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
                Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
                if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                    // First for psi values 
                    if (GlobalBCRS::s_bcLo_psi[i] == 1)
                    {
                        if (!GlobalBCRS::s_printedThatLo_psi[i])
                        {
                            GlobalBCRS::s_printedThatLo_psi[i] = true;
                            pout() << "Constant Neumann bcs imposed on psi for low "
                                      "side direction "
                                   << i << endl;
                        }
                        NeumBC(a_state, valid, a_dx, a_homogeneous,
                               ParseValuePsi, // BCValueHolder class
                               i, Side::Lo, Interval(c_psi, c_psi));
                    }
                    else if (GlobalBCRS::s_bcLo_psi[i] == 0)
                    {
                        if (!GlobalBCRS::s_printedThatLo_psi[i])
                        {
                            GlobalBCRS::s_printedThatLo_psi[i] = true;
                            pout() << "Constant Dirichlet bcs imposed on psi for low "
                                      "side direction "
                                   << i << endl;
                        }
                        DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValuePsi,
                               i, Side::Lo, Interval(c_psi, c_psi));
                    }
                    else
                    {
                        MayDay::Error("bogus bc flag low side psi");
                    }

                    // Now for Vi values 
                    if (GlobalBCRS::s_bcLo_Vi[i] == 1)
                    {
                        if (!GlobalBCRS::s_printedThatLo_Vi[i])
                        {
                            GlobalBCRS::s_printedThatLo_Vi[i] = true;
                            pout() << "Constant Neumann bcs imposed on Vi for low "
                                      "side direction "
                                   << i << endl;
                        }
                        NeumBC(a_state, valid, a_dx, a_homogeneous,
                               ParseValueVi, // BCValueHolder class
                               i, Side::Lo, Interval(c_V1, c_V3));
                    }
                    else if (GlobalBCRS::s_bcLo_Vi[i] == 0)
                    {
                        if (!GlobalBCRS::s_printedThatLo_Vi[i])
                        {
                            GlobalBCRS::s_printedThatLo_Vi[i] = true;
                            pout() << "Constant Dirichlet bcs imposed on Vi for low "
                                      "side direction "
                                   << i << endl;
                        }
                        DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValueVi,
                               i, Side::Lo, Interval(c_V1, c_V3));
                    }
                    else
                    {
                        MayDay::Error("bogus bc flag low side Vi");
                    }
                }

                if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                    // First for psi
                    if (GlobalBCRS::s_bcHi_psi[i] == 1)
                    {
                        if (!GlobalBCRS::s_printedThatHi_psi[i])
                        {
                            GlobalBCRS::s_printedThatHi_psi[i] = true;
                            pout() << "Constant Neumann bcs imposed on psi for high "
                                      "side direction "
                                   << i << endl;
                        }
                        NeumBC(a_state, valid, a_dx, a_homogeneous, ParseValuePsi,
                               i, Side::Hi, Interval(c_psi, c_psi));
                    }
                    else if (GlobalBCRS::s_bcHi_psi[i] == 0)
                    {
                        if (!GlobalBCRS::s_printedThatHi_psi[i])
                        {
                            GlobalBCRS::s_printedThatHi_psi[i] = true;
                            pout() << "Constant Dirichlet bcs imposed on psi for high "
                                      "side direction "
                                   << i << endl;
                        }
                        DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValuePsi,
                               i, Side::Hi, Interval(c_psi, c_psi));
                    }
                    else
                    {
                        MayDay::Error("bogus bc flag high side psi");
                    }

                    // Now for Vi
                    if (GlobalBCRS::s_bcHi_Vi[i] == 1)
                    {
                        if (!GlobalBCRS::s_printedThatHi_Vi[i])
                        {
                            GlobalBCRS::s_printedThatHi_Vi[i] = true;
                            pout() << "Constant Neumann bcs imposed on Vi for high "
                                      "side direction "
                                   << i << endl;
                        }
                        NeumBC(a_state, valid, a_dx, a_homogeneous, ParseValueVi,
                               i, Side::Hi, Interval(c_V1, c_V3));
                    }
                    else if (GlobalBCRS::s_bcHi_Vi[i] == 0)
                    {
                        if (!GlobalBCRS::s_printedThatHi_Vi[i])
                        {
                            GlobalBCRS::s_printedThatHi_Vi[i] = true;
                            pout() << "Constant Dirichlet bcs imposed on psi for high "
                                      "side direction "
                                   << i << endl;
                        }
                        DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValueVi,
                               i, Side::Hi, Interval(c_V1, c_V3));
                    }
                    else
                    {
                        MayDay::Error("bogus bc flag high side Vi");
                    }
                }
            } // else - is periodic

        } // close for idir
    }
}
