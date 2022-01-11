/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef READINUSERVARIABLES_HPP
#define READINUSERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_phi_Re_in, // matter field added
    c_Pi_Re_in,  //(minus) conjugate momentum
    c_phi_Im_in, // matter field added
    c_Pi_Im_in,  //(minus) conjugate momentum

    NUM_READIN_VARS
};

namespace ReadinUserVariables
{
static constexpr char const *variable_names[NUM_READIN_VARS] = {"phi_Re", "Pi_Re", "phi_Im", "Pi_Im"};
}

#endif /* READINUSERVARIABLES_HPP */
