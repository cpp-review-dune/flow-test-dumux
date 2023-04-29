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
 * \brief Exercise to show the spreading of contaminants.
 */

#ifndef DUMUX_LENS2P_PROBLEM_HH
#define DUMUX_LENS2P_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux{

/*!
 * \brief Exercise to show the spreading of contaminants.
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

template <typename TypeTag>
class LensTwoPProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

public:
    /*!
     * \brief The constructor
     */
    LensTwoPProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
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
     *        used for which equation on a given boundary control volume.
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
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
     * \param globalPos the global position
     *
     */
    PrimaryVariables dirichletAtPos( const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        if (onUpperBoundary_(globalPos))
        {
            values[Indices::pressureIdx] = upperPressure_;
            values[Indices::saturationIdx] = 0.0;
        }

        else if (onLowerBoundary_(globalPos))
        {
            values[Indices::pressureIdx] = lowerPressure_;
            values[Indices::saturationIdx] = 0.0;
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     *
     */
    PrimaryVariables neumannAtPos( const GlobalPosition &globalPos) const
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
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     */
    NumEqVector sourceAtPos( const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        // no DNAPL, hydrostatic pressure
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar height = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];

        values[Indices::pressureIdx] = upperPressure_ - depth/height*(upperPressure_-lowerPressure_);
        values[Indices::saturationIdx] = (isInitial_(globalPos)) ?  0.9 :0.0;

        return values;
    }
    // \}

private:
    Scalar eps_ = 3e-6;
    Scalar upperPressure_;
    Scalar lowerPressure_;

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
        return (globalPos[1] > 1.0 + eps_ && globalPos[1] < 1.25 - eps_ && globalPos[0] > 1.7 + eps_
                && globalPos[0] < 2.3 - eps_);
    }
};

} //end namespace Dumux

#endif
