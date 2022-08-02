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
 *
 * \brief The two-phase porousmediumflow problem for exercise-basic
 */

#ifndef DUMUX_EX_BASIC_PROBLEM_2P2C_HH
#define DUMUX_EX_BASIC_PROBLEM_2P2C_HH

#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Gas injection problem where a gas (here  nitrogen) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 *
 * The domain is sized 60 m times 40 m.
 *
 * For the mass conservation equation Neumann boundary conditions are used on
 * the top, on the bottom and on the right of the domain, while Dirichlet conditions
 * apply on the left boundary.
 *
 * Gas is injected at the right boundary from 7 m to 15 m at a rate of
 * 0.001 kg/(s m), the remaining Neumann boundaries are no-flow
 * boundaries.
 *
 * At the Dirichlet boundaries a hydrostatic pressure and a gas saturation of zero a
 *
 * This problem uses the \ref TwoPModel model.
 */
template<class TypeTag>
class Injection2p2cProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    enum { dimWorld = GridView::dimensionworld };
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    Injection2p2cProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                /*tempMax=*/423.15,
                /*numTemp=*/50,
                /*pMin=*/0.0,
                /*pMax=*/30e6,
                /*numP=*/300);

        // name of the problem and output file
        // getParam<TYPE>("GROUPNAME.PARAMNAME") reads and sets parameter PARAMNAME
        // of type TYPE given in the group GROUPNAME from the input file
        name_ = getParam<std::string>("Problem.Name");
        // depth of the aquifer, units: m
        aquiferDepth_ = getParam<Scalar>("Problem.AquiferDepth");
        // the duration of the injection, units: second
        injectionDuration_ = getParam<Scalar>("Problem.InjectionDuration");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    { return name_+"-2p2c"; }

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
         BoundaryTypes bcTypes;
        if (globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        // initialize values to zero, i.e. no-flow Neumann boundary conditions
        NumEqVector values(0.0);

        // if we are inside the injection zone set inflow Neumann boundary conditions
        if (time_ < injectionDuration_
            && globalPos[1] < 15 + eps_ && globalPos[1] > 7 - eps_ && globalPos[0] > 0.9*this->gridGeometry().bBoxMax()[0])
        {
            // TODO: dumux-course-task
            //instead of setting -1e-4 here directly use totalAreaSpecificInflow_ in the computation

            // inject nitrogen. negative values mean injection
            // convert from units kg/(s*m^2) to mole/(s*m^2)
            values[Indices::conti0EqIdx + FluidSystem::N2Idx] = -1e-4/FluidSystem::molarMass(FluidSystem::N2Idx);
            values[Indices::conti0EqIdx + FluidSystem::H2OIdx] = 0.0;
        }

        return values;
    }

    // \}


    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param globalPos The position for which the source term should be evaluated
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(Indices::firstPhaseOnly);
        // get the water density at atmospheric conditions
        const Scalar densityW = FluidSystem::H2O::liquidDensity(this->spatialParams().temperatureAtPos(globalPos), 1.0e5);

        // assume an intially hydrostatic liquid pressure profile
        // note: we subtract rho_w*g*h because g is defined negative
        const Scalar pw = 1.0e5 - densityW*this->spatialParams().gravity(globalPos)[dimWorld-1]*(aquiferDepth_ - globalPos[dimWorld-1]);

        // initially we have some nitrogen dissolved
        // saturation mole fraction would be
        // moleFracLiquidN2 = (pw + pc + p_vap^sat)/henry;
        const Scalar moleFracLiquidN2 = pw*0.95/BinaryCoeff::H2O_N2::henry(this->spatialParams().temperatureAtPos(globalPos));

        // note that because we start with a single phase system the primary variables
        // are pl and x^w_N2. This will switch as soon after we start injecting to a two
        // phase system so the primary variables will be pl and Sn (nonwetting saturation).
        values[Indices::pressureIdx] = pw;
        values[Indices::switchIdx] = moleFracLiquidN2;

        return values;
    }

    // \}

    //! set the time for the time dependent boundary conditions (called from main)
    void setTime(Scalar time)
    { time_ = time; }

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_; //! Problem name
    Scalar aquiferDepth_; //! Depth of the aquifer in m
    Scalar injectionDuration_; //! Duration of the injection in seconds
    Scalar time_;
    //TODO: dumux-course-task
    //define the Scalar totalAreaSpecificInflow_ here

};

} //end namespace Dumux

#endif
