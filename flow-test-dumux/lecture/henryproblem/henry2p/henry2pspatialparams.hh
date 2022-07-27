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
 * \brief The spatial parameters for the HernyProblem which uses the
 *        twophase box model
 */
#ifndef DUMUX_HENRY2P_SPATIAL_PARAMS_HH
#define DUMUX_HENRY2P_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/porousmediumflow/2p/model.hh>

namespace Dumux {

/*!
 * \brief The spatial parameters for the Henry2pProblem which uses the
 *        twophase box model
 */

template<class FVGridGeometry, class Scalar>
class Henry2pSpatialParams :
public FVPorousMediumFlowSpatialParamsMP<FVGridGeometry, Scalar, Henry2pSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = Henry2pSpatialParams<FVGridGeometry, Scalar>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<FVGridGeometry, Scalar, ThisType>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::LinearMaterialDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    Henry2pSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        try
        {
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }
        typename PcKrSwCurve::BasicParams params(0.0/*pcEntry*/, 0.0/*pcMax*/);
        pcKrSwCurve_ = std::make_unique<PcKrSwCurve>(params);

        k_ = 1.019368e-9;
        porosity_=0.35;
        temperature_ = 293.15;
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return k_; }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param globalPos The global position
     * \return the material parameters object
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(*pcKrSwCurve_); }


    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

        /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return temperature_; // in [K]
    };

private:
    PermeabilityType k_;
    Scalar porosity_;
    Scalar temperature_;
    std::unique_ptr<PcKrSwCurve> pcKrSwCurve_;


};

} // end namespace Dumux
#endif
