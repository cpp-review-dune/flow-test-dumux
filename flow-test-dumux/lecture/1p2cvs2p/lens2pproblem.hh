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
 * \brief Soil decontamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENS2P_EXERCISE1_PROBLEM_HH
#define DUMUX_LENS2P_EXERCISE1_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \brief Soil decontamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * DOC ME! The domain is sized 5m times 4m and features a rectangular lens
 * with low permeablility which spans from (0.8 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 50\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 1\,000\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2p -parameterFile test_2p.input</tt>
 */
template <class TypeTag >
class LensTwoPProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    // primary variable indices
    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int saturationIndex = Indices::saturationIdx;
    // equation indices
    static constexpr int contiNEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx;
    // phase indices
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

public:
    /*!
     * \brief The constructor
     *
     */
    LensTwoPProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        //This overwrites the lens settings made in the spatialparameter file
        this->spatialParams().setLensCoords({0.8, 2.0}, {4.0, 3.0});

        // the boundary condition data
        lowerPressure_ = getParam<Scalar>("Boundary.LowerPressure");
        upperPressure_ = getParam<Scalar>("Boundary.UpperPressure");

        infiltrationRate_ = getParam<Scalar>("Boundary.InfiltrationRate");
        infiltrationEndTime_= getParam<Scalar>("Boundary.InfiltrationEndTime");
    }


    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            values.setAllDirichlet();

        else
            values.setAllNeumann();

        if (onInlet_(globalPos))
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The position for which the Dirichlet condition should be evaluated
     *
     */
    PrimaryVariables dirichletAtPos( const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        if (onUpperBoundary_(globalPos))
        {
            values[pressureIdx] = upperPressure_;
            values[saturationIndex] = 0.0;
        }

        else if (onLowerBoundary_(globalPos))
        {
            values[pressureIdx] = lowerPressure_;
            values[saturationIndex] = 0.0;
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     */
    PrimaryVariables neumannAtPos( const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        if (time() <= infiltrationEndTime_ && onInlet_(globalPos))
            values[contiNEqIdx] = -infiltrationRate_; // kg/(m*s)

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
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
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
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos( const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        // no DNAPL, hydrostatic pressure
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar height = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];

        values[pressureIdx] = upperPressure_ - depth/height*(upperPressure_-lowerPressure_);
        values[saturationIndex] = 0.0;

        return values;
    }
    // \}

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
        const Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        const Scalar lambda = (this->gridGeometry().bBoxMax()[0] - globalPos[0])/width;

        return (onUpperBoundary_(globalPos) && (this->gridGeometry().bBoxMax()[0]-0.45*width)/width > lambda + eps_
                                            && (this->gridGeometry().bBoxMax()[0]-0.55*width)/width < lambda - eps_);
    }

    static constexpr Scalar eps_ = 1e-8;
    Scalar episodeLength_;

    Scalar upperPressure_;
    Scalar lowerPressure_;
    Scalar infiltrationRate_;
    Scalar infiltrationEndTime_;
    TimeLoopPtr timeLoop_;
};

} //end namespace Dumux

#endif
