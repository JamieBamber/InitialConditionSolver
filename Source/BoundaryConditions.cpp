/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "RealVect.H"

// Other includes
#include "BoundaryConditions.hpp"
#include <algorithm>
#include <array>
#include <map>
#include <numeric>
#include <string>

// Chombo namespace
#include "UsingNamespace.H"

BoundaryConditions::params_t::params_t()
{
    // set defaults
    hi_boundary.fill(STATIC_BC);
    lo_boundary.fill(STATIC_BC);
    is_periodic.fill(true);
    boundary_conditions_written = false;
    nonperiodic_boundaries_exist = false;
    reflective_boundaries_exist = false;
    vars_parity_multigrid = MultigridUserVariables::vars_parity;
    vars_parity_grchombo = GRChomboUserVariables::vars_parity;
    extrapolation_order = 1;
}

void BoundaryConditions::params_t::set_is_periodic(
    const std::array<bool, CH_SPACEDIM> &a_is_periodic)
{
    is_periodic = a_is_periodic;
    FOR(idir)
    {
        if (!is_periodic[idir])
            nonperiodic_boundaries_exist = true;
    }
}
void BoundaryConditions::params_t::set_hi_boundary(
    const std::array<int, CH_SPACEDIM> &a_hi_boundary)
{
    FOR(idir)
    {
        if (!is_periodic[idir])
        {
            hi_boundary[idir] = a_hi_boundary[idir];
            if (hi_boundary[idir] == REFLECTIVE_BC)
            {
                reflective_boundaries_exist = true;
            }
            else if (hi_boundary[idir] == EXTRAPOLATING_BC)
            {
                extrapolating_boundaries_exist = true;
            }
        }
    }
}
void BoundaryConditions::params_t::set_lo_boundary(
    const std::array<int, CH_SPACEDIM> &a_lo_boundary)
{
    FOR(idir)
    {
        if (!is_periodic[idir])
        {
            lo_boundary[idir] = a_lo_boundary[idir];
            if (lo_boundary[idir] == REFLECTIVE_BC)
            {
                reflective_boundaries_exist = true;
            }
            else if (lo_boundary[idir] == EXTRAPOLATING_BC)
            {
                extrapolating_boundaries_exist = true;
            }
        }
    }
}

void BoundaryConditions::params_t::read_params(GRParmParse &pp)
{
    //    using namespace MultigridUserVariables; // for loading the arrays

    // still load even if not contained, to ensure printout saying parameters
    // were set to their default values
    std::array<bool, CH_SPACEDIM> isPeriodic;
    pp.load("is_periodic", isPeriodic, is_periodic);
    if (pp.contains("is_periodic"))
        set_is_periodic(isPeriodic);

    std::array<int, CH_SPACEDIM> hiBoundary;
    pp.load("hi_boundary", hiBoundary, hi_boundary);
    if (pp.contains("hi_boundary"))
        set_hi_boundary(hiBoundary);

    std::array<int, CH_SPACEDIM> loBoundary;
    pp.load("lo_boundary", loBoundary, lo_boundary);
    if (pp.contains("lo_boundary"))
        set_lo_boundary(loBoundary);

    if (extrapolating_boundaries_exist)
    {
        pp.load("extrapolation_order", extrapolation_order, 1);
    }
    if (nonperiodic_boundaries_exist && !boundary_conditions_written)
    {
        // write out boundary conditions where non periodic - useful for
        // debug
        write_boundary_conditions(*this);
        boundary_conditions_written = true;
    }
}

/// define function sets members and is_defined set to true
void BoundaryConditions::define(double a_dx,
                                std::array<double, CH_SPACEDIM> a_center,
                                const params_t &a_params,
                                ProblemDomain a_domain, int a_num_ghosts)
{
    m_dx = a_dx;
    m_params = a_params;
    m_domain = a_domain;
    m_domain_box = a_domain.domainBox();
    m_num_ghosts = a_num_ghosts;
    FOR(i) { m_center[i] = a_center[i]; }
    is_defined = true;
}

void BoundaryConditions::write_reflective_conditions(int idir,
                                                     const params_t &a_params)
{
    pout() << "The variables that are parity odd in this direction are : "
           << endl;
    for (int icomp = 0; icomp < NUM_MULTIGRID_VARS; icomp++)
    {
        int parity = get_var_parity(icomp, idir, a_params);
        if (parity == -1)
        {
            pout() << MultigridUserVariables::variable_names[icomp] << "    ";
        }
    }
    for (int icomp = 0; icomp < NUM_GRCHOMBO_VARS; icomp++)
    {
        int parity =
            get_var_parity(icomp, idir, a_params, VariableType::grchombo);
        if (parity == -1)
        {
            pout() << GRChomboUserVariables::variable_names[icomp] << "    ";
        }
    }
}

/// write out boundary params (used during setup for debugging)
void BoundaryConditions::write_boundary_conditions(const params_t &a_params)
{
    pout() << "You are using non periodic boundary conditions." << endl;
    pout() << "The boundary params chosen are:  " << endl;
    pout() << "---------------------------------" << endl;

    std::map<int, std::string> bc_names = {{STATIC_BC, "Static"},
                                           {REFLECTIVE_BC, "Reflective"},
                                           {EXTRAPOLATING_BC, "Extrapolating"}};
    FOR(idir)
    {
        if (!a_params.is_periodic[idir])
        {
            pout() << "- " << bc_names[a_params.hi_boundary[idir]]
                   << " boundaries in direction high " << idir << endl;
            // high directions
            if (a_params.hi_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir, a_params);
            }
            pout() << "\n" << endl;

            // low directions
            pout() << "- " << bc_names[a_params.lo_boundary[idir]]
                   << " boundaries in direction low " << idir << endl;
            if (a_params.lo_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir, a_params);
            }
            pout() << "\n" << endl;
        }
    }
    pout() << "---------------------------------" << endl;
}

/// The function which returns the parity of each of the vars in
/// UserVariables.hpp The parity should be defined in the params file, and
/// will be output to the pout files for checking at start/restart of
/// simulation (It is only required for reflective boundary conditions.)
int BoundaryConditions::get_var_parity(int a_comp, int a_dir,
                                       const VariableType var_type) const
{
    int var_parity = get_var_parity(a_comp, a_dir, m_params, var_type);

    return var_parity;
}

/// static version used for initial output of boundary values
int BoundaryConditions::get_var_parity(int a_comp, int a_dir,
                                       const params_t &a_params,
                                       const VariableType var_type)
{
    int comp_parity = (var_type == VariableType::multigrid
                           ? a_params.vars_parity_multigrid[a_comp]
                           : a_params.vars_parity_grchombo[a_comp]);

    int vars_parity = 1;
    if ((a_dir == 0) && (comp_parity == ODD_X || comp_parity == ODD_XY ||
                         comp_parity == ODD_XZ || comp_parity == ODD_XYZ))
    {
        vars_parity = -1;
    }
    else if ((a_dir == 1) && (comp_parity == ODD_Y || comp_parity == ODD_XY ||
                              comp_parity == ODD_YZ || comp_parity == ODD_XYZ))
    {
        vars_parity = -1;
    }
    else if ((a_dir == 2) && (comp_parity == ODD_Z || comp_parity == ODD_XZ ||
                              comp_parity == ODD_YZ || comp_parity == ODD_XYZ))
    {
        vars_parity = -1;
    }
    return vars_parity;
}

/// fill solution boundary conditions, after NL update
void BoundaryConditions::fill_multigrid_boundaries(
    const Side::LoHiSide a_side, LevelData<FArrayBox> &a_state,
    const Interval &a_comps)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions::fill_multigrid_boundaries");

    // cycle through the directions
    FOR(idir)
    {
        // only do something if this direction is not periodic and solution
        // boundary enforced in this direction
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);

            fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                    boundary_condition, a_comps,
                                    VariableType::multigrid);
        }
    }
}

/// fill grchombo boundary conditions, after solution found
void BoundaryConditions::fill_grchombo_boundaries(const Side::LoHiSide a_side,
                                                  LevelData<FArrayBox> &a_state,
                                                  const Interval &a_comps)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions::fill_grchombo_boundaries");

    // cycle through the directions
    FOR(idir)
    {
        // only do something if this direction is not periodic and solution
        // boundary enforced in this direction
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);

            fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                    boundary_condition, a_comps,
                                    VariableType::grchombo);
        }
    }
}

/// Fill the boundary values appropriately based on the params set
/// in the direction dir
void BoundaryConditions::fill_boundary_cells_dir(
    const Side::LoHiSide a_side, const LevelData<FArrayBox> &a_soln,
    LevelData<FArrayBox> &a_out, const int dir, const int boundary_condition,
    const Interval &a_comps, const VariableType var_type)
{
    std::vector<int> comps_vector;
    comps_vector.resize(a_comps.size());
    std::iota(comps_vector.begin(), comps_vector.end(), a_comps.begin());

    // iterate through the boxes, shared amongst threads
    DataIterator dit = a_out.dataIterator();
    int nbox = dit.size();
#pragma omp parallel for default(shared)
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex dind = dit[ibox];
        FArrayBox &out_box = a_out[dind];
        const FArrayBox &soln_box = a_soln[dind];
        Box this_box = out_box.box();
        IntVect offset_lo = -this_box.smallEnd() + m_domain_box.smallEnd();
        IntVect offset_hi = +this_box.bigEnd() - m_domain_box.bigEnd();

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        this_box &= m_domain_box;
        // get the boundary box (may be Empty)
        Box boundary_box =
            get_boundary_box(a_side, dir, offset_lo, offset_hi, this_box);

        // now we have the appropriate box, fill it!
        BoxIterator bit(boundary_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            switch (boundary_condition)
            {
            // simplest case - boundary values are set to zero
            case STATIC_BC:
            {
                for (int icomp = a_comps.begin(); icomp <= a_comps.end();
                     ++icomp)
                {
                    out_box(iv, icomp) = 0.0;
                }
                break;
            }
            // Enforce a reflective symmetry in some direction
            case REFLECTIVE_BC:
            {
                fill_reflective_cell(out_box, iv, a_side, dir, comps_vector,
                                     var_type);
                break;
            }
            case EXTRAPOLATING_BC:
            {
                fill_extrapolating_cell(out_box, iv, a_side, dir, comps_vector,
                                        m_params.extrapolation_order);
                break;
            }
            default:
                MayDay::Error(
                    "BoundaryCondition::Supplied boundary not supported.");
            } // end switch
        }     // end iterate over box
    }         // end iterate over boxes
}

void BoundaryConditions::fill_reflective_cell(
    FArrayBox &out_box, const IntVect iv, const Side::LoHiSide a_side,
    const int dir, const std::vector<int> &reflective_comps,
    const VariableType var_type) const
{
    // assume boundary is a reflection of values within the grid
    // care must be taken with variable parity to maintain correct
    // values on reflection, e.g. x components of vectors are odd
    // parity in the x direction
    IntVect iv_copy = iv;
    /// where to copy the data from - mirror image in domain
    if (a_side == Side::Lo)
    {
        iv_copy[dir] = -iv[dir] - 1;
    }
    else
    {
        iv_copy[dir] = 2 * m_domain_box.bigEnd(dir) - iv[dir] + 1;
    }

    // replace value at iv with value at iv_copy
    for (int icomp : reflective_comps)
    {
        int parity = get_var_parity(icomp, dir, var_type);
        out_box(iv, icomp) = parity * out_box(iv_copy, icomp);
    }
}

void BoundaryConditions::fill_extrapolating_cell(
    FArrayBox &out_box, const IntVect iv, const Side::LoHiSide a_side,
    const int dir, const std::vector<int> &extrapolating_comps,
    const int order) const
{
    for (int icomp : extrapolating_comps)
    {
        // current radius
        double radius =
            get_radius(iv, m_dx, {m_center[0], m_center[1], m_center[2]});

        // vector of 2 nearest values and radii within the grid
        std::array<double, 2> value_at_point;
        std::array<double, 2> r_at_point;
        // how many units are we from domain boundary?
        int units_from_edge = 0;
        if (a_side == Side::Hi)
        {
            // how many units are we from domain boundary?
            units_from_edge = iv[dir] - m_domain_box.bigEnd(dir);
            // vector of 2 nearest values and radii within the grid
            for (int i = 0; i < 2; i++)
            {
                IntVect iv_tmp = iv;
                iv_tmp[dir] += -units_from_edge - i;
                FOR(idir)
                {
                    if (iv_tmp[idir] > m_domain_box.bigEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.bigEnd(idir);
                    }
                    else if (iv_tmp[idir] < m_domain_box.smallEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.smallEnd(idir);
                    }
                }
                value_at_point[i] = out_box(iv_tmp, icomp);
                r_at_point[i] = get_radius(
                    iv_tmp, m_dx, {m_center[0], m_center[1], m_center[2]});
            }
        }
        else // Lo side
        {
            // how many units are we from domain boundary?
            units_from_edge = -iv[dir] + m_domain_box.smallEnd(dir);
            // vector of 2 nearest values within the grid
            for (int i = 0; i < 2; i++)
            {
                IntVect iv_tmp = iv;
                iv_tmp[dir] += units_from_edge + i;
                FOR(idir)
                {
                    if (iv_tmp[idir] > m_domain_box.bigEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.bigEnd(idir);
                    }
                    else if (iv_tmp[idir] < m_domain_box.smallEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.smallEnd(idir);
                    }
                }
                value_at_point[i] = out_box(iv_tmp, icomp);
                r_at_point[i] = get_radius(
                    iv_tmp, m_dx, {m_center[0], m_center[1], m_center[2]});
            }
        }

        // assume some radial dependence and fit it
        double analytic_change = 0.0;
        // comp = const
        if (order == 0)
        {
            analytic_change = 0.0;
        }
        // comp = B + A*r
        else if (order == 1)
        {
            double delta_r_in_domain = r_at_point[1] - r_at_point[0];
            double A =
                (value_at_point[1] - value_at_point[0]) / delta_r_in_domain;
            double delta_r_here = radius - r_at_point[0];
            analytic_change = A * delta_r_here;
        }
        // other orders not supported yet
        else
        {
            MayDay::Error("Order not supported for boundary extrapolation.");
        }

        // set the value here to the extrapolated value
        out_box(iv, icomp) = value_at_point[0] + analytic_change;
    }
}

/// Get the boundary condition for a_dir and a_side
int BoundaryConditions::get_boundary_condition(const Side::LoHiSide a_side,
                                               const int a_dir)
{
    int boundary_condition = 0;
    if (a_side == Side::Lo)
    {
        boundary_condition = m_params.lo_boundary[a_dir];
    }
    else
    {
        boundary_condition = m_params.hi_boundary[a_dir];
    }
    return boundary_condition;
}

/// get the boundary box to fill if we are at a boundary
Box BoundaryConditions::get_boundary_box(
    const Side::LoHiSide a_side, const int a_dir, const IntVect &offset_lo,
    const IntVect &offset_hi, Box &this_ghostless_box, int shrink_for_coarse)
{
    // default constructor gives empty box
    Box boundary_box;

    // check if we are over the edges of the domain - are we a boundary box?
    // if so create the box of the cells we want to fill
    if (((a_side == Side::Hi) && (offset_hi[a_dir] > 0)) ||
        ((a_side == Side::Lo) && (offset_lo[a_dir] > 0)))
    {
        // Get just the boundary box to iterate over, m_num_ghosts ghost
        // cells unless we are filling the coarse cells in the interp case
        // where we want to fill only two coarse ghost cells (to cover 3
        // fine ones)
        if (a_side == Side::Lo)
        {
            boundary_box = adjCellLo(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }
        else
        {
            boundary_box = adjCellHi(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }

        // adjust for any offsets - catches the corners etc
        // but only want to fill them once, so y fills x, z fills y and x
        // etc. Required in periodic direction corners in cases where there
        // are mixed boundaries, (otherwise these corners are full of nans)
        FOR(idir)
        {
            if (offset_lo[idir] > 0) // this direction is a low end boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growLo(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growLo(idir, -shrink_for_coarse);
                }
            }

            if (offset_hi[idir] > 0) // this direction is a high end
                                     // boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growHi(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growHi(idir, -shrink_for_coarse);
                }
            }
        }
    }
    return boundary_box;
}

/// Operator called by transform to grow the boxes where required
Box ExpandGridsToBoundaries::operator()(const Box &a_in_box)
{
    Box out_box = a_in_box;
    IntVect offset_lo =
        -out_box.smallEnd() + m_boundaries.m_domain_box.smallEnd();
    IntVect offset_hi = +out_box.bigEnd() - m_boundaries.m_domain_box.bigEnd();

    FOR(idir)
    {
        if (!m_boundaries.m_params.is_periodic[idir])
        {
            if (offset_lo[idir] == 0)
            {
                out_box.growLo(idir, m_boundaries.m_num_ghosts);
            }
            if (offset_hi[idir] == 0)
            {
                out_box.growHi(idir, m_boundaries.m_num_ghosts);
            }
        }
    }
    return out_box;
}

/// This function takes a default constructed open DisjointBoxLayout and
/// grows the boxes lying along the boundary to include the boundaries if
/// necessary.
/// It is used to define the correct
/// DisjointBoxLayout for the exchange copier so that shared boundary ghosts are
/// exchanged correctly.
void BoundaryConditions::expand_grids_to_boundaries(
    DisjointBoxLayout &a_out_grids, const DisjointBoxLayout &a_in_grids)
{
    if (!a_in_grids.isClosed())
    {
        MayDay::Error("input to expand_grids_to_boundaries must be closed");
    }

    // Grow the problem domain to include the boundary ghosts
    ProblemDomain domain_with_boundaries = a_in_grids.physDomain();
    FOR(idir)
    {
        if (!m_params.is_periodic[idir])
        {
            domain_with_boundaries.growLo(idir, m_num_ghosts);
            domain_with_boundaries.growHi(idir, m_num_ghosts);
        }
    }

    // Copy grids and apply transformation
    a_out_grids.deepCopy(a_in_grids, domain_with_boundaries);
    ExpandGridsToBoundaries expand_grids_to_boundaries(*this);
    a_out_grids.transform(expand_grids_to_boundaries);
    a_out_grids.close();
}

double
BoundaryConditions::get_radius(IntVect integer_coords, double dx,
                               std::array<double, CH_SPACEDIM> center) const
{
    double xx = (integer_coords[0] + 0.5) * dx - center[0];
    double yy = (integer_coords[1] + 0.5) * dx - center[1];
    double zz = (integer_coords[2] + 0.5) * dx - center[2];

    double r = sqrt(xx * xx + yy * yy + zz * zz);

    const double minimum_r = 1e-6;
    return max(r, minimum_r);
}
