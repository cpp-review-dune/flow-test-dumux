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
 * \brief Exercise to show the diffusive spreading of contaminants.
 */

#ifndef DUMUX_LENS_1P2C_PROBLEM_HH
#define DUMUX_LENS_1P2C_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \brief Exercise to show the diffusive spreading of contaminants.
 *
 * The whole problem is symmetric.
 *
 * The domain is sized 2m times 4m and features a rectangular lens
 * with low permeablility which is located in the lower half of the domain.
  *
 * On the top and the bottom of the domain Dirichlet boundary conditions
 * are used, while Neumann conditions apply on the left and right
 * boundaries.
 *
 * The contaminant is initially located in the middle of the domain,
 * but only in the upper layer.
 *
 * The Dirichlet boundaries on the top boundary can be chosen
 * differently to the bottom pressure.
 *
 * Typical simulation parameters can be found in the input file:
 * exercise3.input
 *
 */
template <class TypeTag >
class LensOnePTwoCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    enum {
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx)
    };
public:
    /*!
     * \brief The constructor
     * \param fvGridGeometry The finite volume grid geometry
     */
    LensOnePTwoCProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry )
    : ParentType(fvGridGeometry)
    {
        // the boundary condition data
        lowerPressure_ = getParam<Scalar>("Boundary.LowerPressure");
        upperPressure_ = getParam<Scalar>("Boundary.UpperPressure");

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/ 273.15 + getParam<Scalar>("SpatialParams.Temperature") - 1.0,
        /*Tmax=*/ 273.15 + getParam<Scalar>("SpatialParams.Temperature") + 1.0,
        /*nT=*/3,
        /*pmin=*/1e5,
        /*pmax=*/3e7,
        /*np=*/200);
    }


    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            values.setAllDirichlet();

        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The position for which the Dirichlet condition should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
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
     *
     * \param globalPos The position for which the neumann condition should be evaluated
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
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
     * \param globalPos The position for which the source should be evaluated
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     */
    PrimaryVariables initialAtPos( const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        // no contaminant, hydrostatic pressure
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar height = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];

        values[Indices::conti0EqIdx] = upperPressure_ - depth/height*(upperPressure_ - lowerPressure_);
        values[N2Idx] = (isInitial_(globalPos)) ? 9.12825137e-6 : 0.0;

        return values;
    }
    // \}

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
    bool isInitial_(const GlobalPosition &globalPos) const
    {
        return (globalPos[1] > 1.0 && globalPos[1] < 1.25 && globalPos[0] > 1.7
                && globalPos[0] < 2.3);
    }

    static constexpr Scalar eps_ = 1e-8;

    Scalar episodeLength_;

    Scalar upperPressure_;
    Scalar lowerPressure_;
    Scalar infiltrationRate_;
    Scalar infiltrationStartTime_;
    Scalar infiltrationEndTime_;
    bool paraviewOutput_;
    std::vector<int> numCells_;

};

} //end namespace Dumux

#endif
