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
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */

#ifndef DUMUX_EXGRIDS_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_EXGRIDS_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotpckrsw.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */
template<class GridGeometry, class Scalar>
class InjectionSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, InjectionSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = InjectionSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;

    // get the dimensions of the simulation domain from GridView
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param gridGeometry The finite volume grid geometry
     */
    InjectionSpatialParams(std::shared_ptr<const GridGeometry>& gridGeometry)
    : ParentType(gridGeometry)
    , aquitardPcKrSwCurve_("SpatialParams.Aquitard")
    , aquiferPcKrSwCurve_("SpatialParams.Aquifer")
    {
        // Aquifer Height, measured from the bottom
        aquiferHeightFromBottom_ = 30.0;

        // intrinsic permeabilities
        aquitardK_ = getParam<Scalar>("SpatialParams.PermeabilityAquitard");
        aquiferK_ = getParam<Scalar>("SpatialParams.PermeabilityAquifer");

        // porosities
        aquitardPorosity_ = 0.2;
        aquiferPorosity_ = 0.4;
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        // here, either aquitard or aquifer permeability are returned, depending on the global position
        if (isInAquitard_(globalPos))
            return aquitardK_;
        return aquiferK_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        // here, either aquitard or aquifer porosity are returned, depending on the global position
        if (isInAquitard_(globalPos))
            return aquitardPorosity_;
        return aquiferPorosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isInAquitard_(globalPos))
            return makeFluidMatrixInteraction(aquitardPcKrSwCurve_);
        return makeFluidMatrixInteraction(aquiferPcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The position of the center of the element
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

    /*!
     * \brief Returns the temperature at the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 30;
    }

private:

    static constexpr Scalar eps_ = 1e-6;

    // provides a convenient way distinguishing whether a given location is inside the aquitard
    bool isInAquitard_(const GlobalPosition &globalPos) const
    {
        // globalPos[dimWorld-1] is the y direction for 2D grids or the z direction for 3D grids
        return globalPos[dimWorld-1] > aquiferHeightFromBottom_ + eps_;
    }

    Scalar aquitardK_;
    Scalar aquiferK_;
    Scalar aquiferHeightFromBottom_;


    Scalar aquitardPorosity_;
    Scalar aquiferPorosity_;

    const PcKrSwCurve aquitardPcKrSwCurve_;
    const PcKrSwCurve aquiferPcKrSwCurve_;
};

} // end namespace Dumux

#endif
