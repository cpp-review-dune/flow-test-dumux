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
 * \brief Definition of the spatial parameters for the lens
 *        problem which uses the isothermal 2p or 2p2c box model
 */

#ifndef DUMUX_LENS2P_SPATIALPARAMS_HH
#define DUMUX_LENS2P_SPATIALPARAMS_HH

#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>

namespace Dumux {


/*!
 * \brief Definition of the spatial parameters for the lens
 *        problem which uses the isothermal 2p or 2p2c box model
 */
template<class FVGridGeometry, class Scalar>
class Lens2pSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<FVGridGeometry, Scalar, Lens2pSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ThisType = Lens2pSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<FVGridGeometry, Scalar, ThisType>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSw = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     */
    Lens2pSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    , pcKrSwFine_("SpatialParams.Fine")
    , pcKrSwCoarse_("SpatialParams.Coarse")
    {
        lensLowerLeft_ = {0.0, 0.0};
        lensUpperRight_= {4.0, 1.0};

        lensPorosity_ = getParam<Scalar> ("SpatialParams.Fine.Porosity");
        outerPorosity_ = getParam<Scalar>("SpatialParams.Coarse.Porosity");

        lensK_ = getParam<Scalar>("SpatialParams.Fine.Permeability");
        outerK_ = getParam<Scalar>("SpatialParams.Coarse.Permeability");

        temperature_ = getParam<Scalar>("SpatialParams.Temperature");
    }
    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
            return lensK_;
         return outerK_;
    }

    /*!
     * \brief Function for defining the porosity.
     *        That is possibly solution dependent.
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
            return lensPorosity_;
        return outerPorosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
            return makeFluidMatrixInteraction(pcKrSwFine_);
        return makeFluidMatrixInteraction(pcKrSwCoarse_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    /*!
     * \brief Returns the temperature at the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return temperature_;
    }

    /*!
     * \brief Set the bounding box of the fine-sand lens
     * \param lensLowerLeft the lower left corner coordinates vector
     * \param lensUpperRight the upper right corner coordinates vector
     */
    void setLensCoords(const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
    { lensLowerLeft_ = lowerLeft; lensUpperRight_ = upperRight;}

private:
    /*!
     * \brief Return if a point is inside the lens domain
     * \param globalPos The global position of the point to check
     */
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        // points on the lens boundary are considered inside
        for (int i = 0; i < dimWorld; ++i)
            if (globalPos[i] < lensLowerLeft_[i] - eps_ || globalPos[i] > lensUpperRight_[i] + eps_)
                return false;
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;
    Scalar lensPorosity_;
    Scalar outerPorosity_;
    Scalar temperature_;
    const PcKrSw pcKrSwFine_;
    const PcKrSw pcKrSwCoarse_;
    const Scalar eps_{1e-6};
};

} // end namespace Dumux

#endif
