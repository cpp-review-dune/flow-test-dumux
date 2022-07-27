// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
#ifndef DUMUX_NAPLINFILTRATIONPROBLEM_3P_3C_HH
#define DUMUX_NAPLINFILTRATIONPROBLEM_3P_3C_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux
{

/*!
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep, where
 * the lower part represents a slightly inclined groundwater table, and the
 * upper part is the vadose zone.
 * A LNAPL (Non-Aqueous Phase Liquid which is lighter than water) infiltrates
 * (modelled with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and migrates
 * on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually trapped
 * immobile NAPL, which can in the following dissolve and evaporate slowly,
 * and eventually be transported by advection and diffusion.
 *
 * Left and right boundaries are constant head boundaries (Dirichlet),
 * Top and bottom are Neumann boundaries, all no-flow except for the small
 * infiltration zone in the upper left part.
 *
 * This problem uses the \ref ThreePThreeCModel.
 *
 * This problem should typically be simulated for 30 days.
 * A good choice for the initial time step size is 60 s.
 * To adjust the simulation time it is necessary to edit the file test_box3p3c.input
 * or test_cc3p3c.input.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box3p3c test_box3p3c.input</tt> or
 * <tt>./test_cc3p3c test_cc3p3c.input</tt>
 *  */
template <class TypeTag >
class InfiltrationThreePThreeCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    // copy some indices for convenience
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        // phase state
        wgPhaseOnly = Indices::wgPhaseOnly,

        //!< Index of the mass conservation equation for the contaminant component
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nCompIdx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InfiltrationThreePThreeCProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        FluidSystem::init(282.15, 284.15, 3, 8e4, 3e5, 200);

        name_ = getParam<std::string>("Problem.Name");
        time_ = 0.0;
    }

    void setTime(Scalar time)
    { time_ = time; }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_
           || globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        // negative values for injection
        if (time_ < 2592000.0 - eps_)
        {
            if ((globalPos[0] < 175.0 + eps_) && (globalPos[0] > 150.0 - eps_)
                 && (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_))
            {
                //mole flow conversion to mass flow with molar mass M(Mesit.)=0,120 kg/mol --> 1.2e-4 kg/(sm)
                //the 3p3c model uses mole fractions
                values[contiNEqIdx] = -0.001;
            }
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

private:
    // internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(wgPhaseOnly);
        Scalar swr=0.12, sgr=0.03;

        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        if (y > (-1e-3*x + 5) - eps_)
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5e-4*x + 5));
            if (pc < 0.0) pc = 0.0;

            Scalar sw = invertPcgw_(pc, this->spatialParams().fluidMatrixInteractionAtPos(globalPos));
            if (sw < swr)
                sw = swr;
            if (sw > 1.-sgr)
                sw = 1.-sgr;

            values[pressureIdx] = 1e5 ;
            values[switch1Idx] = sw;
            values[switch2Idx] = 0.0;
        }
        else
        {
            values[pressureIdx] = 1e5 + 9.81 * 1000.0 * ((-5e-4*x + 5) - y);
            values[switch1Idx] = 1.-sgr;
            values[switch2Idx] = 0.0;
        }
        return values;
    }

    // small solver inverting the pc curve
    template<class FMInteraction>
    static Scalar invertPcgw_(Scalar pcIn, const FMInteraction& fluidMatrixInteraction)
    {
        Scalar lower(0.0);
        Scalar upper(1.0);
        const unsigned int maxIter = 25;
        const Scalar bisLimit = 1.0;

        Scalar sw, pcgw;
        for (unsigned int k = 1; k <= maxIter; k++)
        {
            sw = 0.5*(upper + lower);
            pcgw = fluidMatrixInteraction.pcgw(sw, 0.0/*sn*/);
            const Scalar delta = std::abs(pcgw - pcIn);
            if (delta < bisLimit)
                return sw;

            if (k == maxIter)
                return sw;

            if (pcgw > pcIn)
                lower = sw;
            else
                upper = sw;
        }
        return sw;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    Scalar time_;
};

} // end namespace Dumux

#endif
