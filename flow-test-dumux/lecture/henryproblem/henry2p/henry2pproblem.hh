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
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_HENRY2PPROBLEM_HH
#define DUMUX_HENRY2PPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template <class TypeTag >
class Henry2pProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    // primary variable indices
    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int saturationIdx = Indices::saturationIdx;
    // equation indices
    static constexpr int conti0EqIdx = Indices::conti0EqIdx;
    // Grid and world dimension
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    /*!
     * \brief The constructor
     *
     */
    Henry2pProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        FluidSystem::init();
        eps_ = 3e-6;
        freshWaterFluxRate_= getParam<Scalar>("Problem.freshWaterFluxRate");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    NumEqVector sourceAtPos( const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(onRightBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        if (globalPos[0]>2-eps_)// if(onRightBoundary_(globalPos))
        {

            if(globalPos[1]<0.8+eps_)
            {
                Scalar densityB = 1025;
                values[pressureIdx] = 1.0133e5+(depthBOR_-globalPos[1]-0.2)*densityB*9.81+1000*9.81*0.2;
                values[saturationIdx] = 1.0;
            }

            else
            {
                Scalar densityB = 1025;
                values[pressureIdx] = 1.0133e5+(depthBOR_-globalPos[1]-0.2)*densityB*9.81+1000*9.81*0.2;
                values[saturationIdx] = 0.0;
            }
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        values = 0.0;
        if (onLeftBoundary_(globalPos))
        {
            values[conti0EqIdx] = -freshWaterFluxRate_;//-6.6E-2;// kg / (m * s)
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
     * \param values The initial values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos( const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        if(globalPos[1]<0.8+eps_)
        {
            Scalar densityB = 1025;
            values[pressureIdx] = 1.0133e5+(depthBOR_-globalPos[1]-0.2)*densityB*9.81+1000*9.81*0.2;
            values[saturationIdx] = 0.0;
        }

        else
        {
            Scalar densityW = 1000;
            values[pressureIdx] = 1.0133e5+(depthBOR_-globalPos[1])*densityW*9.81;
            values[saturationIdx] = 0.0;
        }

        return values;
    }
    // \}

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    static constexpr  Scalar depthBOR_ = 1.0;
    Scalar eps_;
    Scalar freshWaterFluxRate_;
};

} //end namespace Dumux

#endif
