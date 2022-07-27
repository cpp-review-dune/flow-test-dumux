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
/**
 * \file
 * \brief Definition of the 1p2c Henry problem
  *
 */
#ifndef DUMUX_HENRY1P2C_PROBLEM_HH
#define DUMUX_HENRY1P2C_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template <class TypeTag>
class Henry1p2cProblem : public PorousMediumFlowProblem<TypeTag>
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

     // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
      // indices of the equations
    static constexpr int conti0EqIdx = Indices::conti0EqIdx;

    // Grid and world dimension
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    Henry1p2cProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
      FluidSystem::init();
      freshWaterFluxRate_= getParam<Scalar>("Problem.freshWaterFluxRate");
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();

        if (globalPos[0]<eps_  || globalPos[1]<eps_ || globalPos[1]>1-eps_)
        {
            values.setAllNeumann();
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        initial_(values, globalPos);
        if (globalPos[0]>2-eps_ && globalPos[1]<0.8+eps_)
        {
            values[FluidSystem::SaltIdx] =  0.03922;
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        const auto globalPos = element.geometry().corner(scvf.insideScvIdx());

        NumEqVector values(0.0);
        if (globalPos[0]<eps_)
        {
            values[conti0EqIdx] = -freshWaterFluxRate_; //-6.6E-2; //6.6e-5*m2/s*1000kg/m^3
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
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        initial_(values, globalPos);

        return values;
    }
    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1025;
        values[Indices::pressureIdx] = 1.0133e5+(depthBOR_-globalPos[1])*densityW*9.81; //initial condition for the pressure
        values[FluidSystem::SaltIdx] = 0.0; //initial condition for the salt molefraction
    }

    Scalar freshWaterFluxRate_;
    static constexpr Scalar eps_ = 1e-6;
    static constexpr Scalar depthBOR_ = 1.0; // [m]
};

} // end namespace Dumux

#endif
