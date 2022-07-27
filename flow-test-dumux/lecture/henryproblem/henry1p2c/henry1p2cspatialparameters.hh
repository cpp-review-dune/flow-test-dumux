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
 * \brief Definition of the spatial parameters for the 1p2c
 *       henry problem
 */
#ifndef DUMUX_HENRY1P2C_SPATIAL_PARAMETERS_HH
#define DUMUX_HENRY1P2C_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 *
 * \brief Definition of the spatial parameters for the 1p2c
 *        Henry problem
 */
template<class FVGridGeometry, class Scalar>
class Henry1p2cSpatialParams :
public FVPorousMediumFlowSpatialParamsOneP<FVGridGeometry, Scalar,
                       Henry1p2cSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<FVGridGeometry, Scalar,
                                       Henry1p2cSpatialParams<FVGridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    using PermeabilityType = Scalar;

    Henry1p2cSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // intrinsic permeabilities
        perm_ = 1.019368e-9;
        porosity_ = 0.35;
        temperature_ = 293.15;
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return perm_; }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    double porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Defines the dispersion tensor \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    std::array<Scalar, 2> dispersionAlphas(const GlobalPosition& globalPos,
                                           int phaseIdx = 0,
                                           int compIdx = 0) const
    {
        static const auto alphaL = getParam<Scalar>("Problem.AlphaL");
        static const auto alphaT = getParam<Scalar>("Problem.AlphaT");
        return std::array<Scalar, 2> {alphaL, alphaT};
    }

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
    PermeabilityType perm_;
    Scalar porosity_;
    Scalar temperature_;
};

} // end namespace Dumux

#endif
