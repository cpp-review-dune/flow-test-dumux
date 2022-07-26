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
 * \brief Soil decontamination problem where a contaminant infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENS_1P2C_PROBLEM_HH
#define DUMUX_LENS_1P2C_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \brief Soil decontamination problem where a contaminant infiltrates a fully
 *        water saturated medium.
 *
 * The whole problem is symmetric.
 *
 * The domain is sized 5m times 4m and features a rectangular lens
 * with low permeablility which spans from (0.8 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability.
 *
 * On the top and the bottom of the domain Dirichlet boundary conditions
 * are used, while Neumann conditions apply on the left and right
 * boundaries.
 *
 * The contaminant is injected at the top boundary from 2.25m to 2.75m at a variable rate.
 *
 * The Dirichlet boundaries on the top boundary is different to the bottom pressure.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically simulated until \f$t_{\text{end}} =
 * 50\,000\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 1\,000\;s\f$.
 *
 */
template <class TypeTag>
class LensOnePTwoCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;
    enum {
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx)
    };
public:
    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The Finite-Volume-Grid-Geometry
     */
    LensOnePTwoCProblem( std::shared_ptr<const FVGridGeometry> fvGridGeometry )
    : ParentType(fvGridGeometry)
    {
        FluidSystem::init();

        //This overwrites the lens settings made in the spatialparameter file
        this->spatialParams().setLensCoords({0.8, 2.0}, {4.0, 3.0});

        infiltrationRate_ = getParam<Scalar>("Boundary.InfiltrationRate");
        infiltrationStartTime_= 1.0e-9; //The infiltrations starts always after the first time step!
        infiltrationEndTime_= getParam<Scalar>("Boundary.InfiltrationEndTime");

        // the boundary condition data
        lowerPressure_ = getParam<Scalar>("Boundary.LowerPressure");
        upperPressure_ = getParam<Scalar>("Boundary.UpperPressure");
    }



    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos DOC ME!
     *
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            values.setAllDirichlet();

        else
            values.setAllNeumann();

        if (onInlet_(globalPos))
            values.setNeumann(N2Idx);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     * \param globalPos he position for which the dirichlet condition should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos( const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        if (onUpperBoundary_(globalPos))
        {
            values[Indices::pressureIdx] = upperPressure_;
        }

        else if (onLowerBoundary_(globalPos))
        {
            values[Indices::pressureIdx] = lowerPressure_;
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     * \param globalPos he position for which the neumann condition should be evaluated
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        if (time() <= infiltrationEndTime_ && infiltrationStartTime_ <= time() && onInlet_(globalPos))
            values[N2Idx] = -infiltrationRate_;

        return values;
    }

    /*!
     * \name Volume terms
     */

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     * \param globalPos he position for which the source should be evaluated
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos the position for which the initial condition should be evaluated
     *
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        // no contaminant, hydrostatic pressure
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar height = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];

        values[Indices::pressureIdx] = upperPressure_ - depth/height*(upperPressure_-lowerPressure_);

        return values;
    }

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

private:
    /*!
     * \brief Returns true if the point is located on the lower boundary
     * \param globalPos The position of the point in global coordinates
     */
    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }
    /*!
     * \brief Returns true if the point is located on the upper boundary
     * \param globalPos The position of the point in global coordinates
     */
    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    /*!
     * \brief returns true if point is located between the given points
     * \param globalPos The position of the point in global coordinates
     */
    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        Scalar lambda = (this->gridGeometry().bBoxMax()[0] - globalPos[0])/width;

        return (onUpperBoundary_(globalPos) && (this->gridGeometry().bBoxMax()[0]-0.45*width)/width > lambda
                                            && lambda > (this->gridGeometry().bBoxMax()[0]-0.55*width)/width);
    }

    Scalar eps_ = 1e-8;
    Scalar episodeLength_;

    Scalar upperPressure_;
    Scalar lowerPressure_;
    Scalar infiltrationRate_;
    Scalar infiltrationStartTime_;
    Scalar infiltrationEndTime_;
    bool paraviewOutput_;
    TimeLoopPtr timeLoop_;
    std::vector<int> numCells_;
};

} //end namespace Dumux

#endif
