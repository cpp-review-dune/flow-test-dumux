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
 * \brief Definition of the spatial parameters for the heterogeneous
 *        problem which uses the non-isothermal or isothermal CO2
 *        fully implicit model.
 */

#ifndef DUMUX_PLUMESHAPE_SPATIAL_PARAMS_HH
#define DUMUX_PLUMESHAPE_SPATIAL_PARAMS_HH

#include <dumux/io/grid/griddata.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

#include <dumux/porousmediumflow/co2/model.hh>

namespace Dumux {

/*!
 * \brief Definition of the spatial parameters for the CO2 plumeshape
 *        problem which uses the non-isothermal or isothermal CO2
 *        fully implicit model.
 */
template<class FVGridGeometry, class Scalar>
class PlumeShapeSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<FVGridGeometry,
                         Scalar,
                         PlumeShapeSpatialParams<FVGridGeometry, Scalar>>
{
    using Grid = typename FVGridGeometry::Grid;
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<FVGridGeometry, Scalar,
                                       PlumeShapeSpatialParams<FVGridGeometry, Scalar>>;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     */
    PlumeShapeSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {
        // intrinsic permeability
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");

        // porosity
        porosity_ = getParam<Scalar>("SpatialParams.Porosity");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return instrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    { return permeability_; }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element
     * \return porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    { return porosity_; }


    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \return the material parameters object
     * \param globalPos The position of the center of the element
     */
    const auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return pcKrSwCurve_; }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The position of the center of the element
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::BrineIdx; }

private:

    Scalar permeability_;
    Scalar porosity_;
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
