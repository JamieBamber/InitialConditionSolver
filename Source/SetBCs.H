#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SETBCS_H_
#define _SETBCS_H_

#include <iostream>
using std::cerr;

#include "AMRMultiGrid.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "DimensionDefinitions.hpp"
#include "VariableCoeffPoissonOperator.H"
#include "MultigridUserVariables.hpp"

// BC class for the linear solver steps
class GlobalBCRS
{
  public:
    static std::vector<int> s_bcLo, s_bcHi;
    static bool s_areBCsParsed, s_valueParsed_single, s_valueParsed_vector;
};

// Allows read in of values to apply in psi BCs
extern void ParseSingleValue(Real *pos, int *dir, Side::LoHiSide *side,
                       Real *a_values);

// Allows read in of values to apply in Vi BCs
extern void ParseVectorValue(Real *pos, int *dir, Side::LoHiSide *side,
                       Real *a_values);

// This is the function that the solver looks for to apply the BCs to the 
// variables in a_state
extern void ParseBC(FArrayBox &a_state, const Box &a_valid,
                    const ProblemDomain &a_domain, Real a_dx,
                    bool a_homogeneous);

extern void ZeroBC(FArrayBox&           a_state,
                   const Box&           a_valid,
                   const ProblemDomain& a_domain,
                   Real                 a_dx,
                   bool                 a_homogeneous,
                   Interval             a_interval);

#endif
