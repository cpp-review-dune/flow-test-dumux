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
 * \brief Properties file of CO2 injection into a reservoir.
 */
#ifndef DUMUX_PLUMESHAPE_PROPERTIES_HH
#define DUMUX_PLUMESHAPE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/co2/model.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/brineco2.hh>

#include <lecture/common/co2tablesbenchmark3.hh>

#include "co2plumeshapespatialparameters.hh"
#include "co2plumeshapeproblem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PlumeShapeTypeTag { using InheritsFrom = std::tuple<TwoPTwoCCO2NI>; };
struct PlumeShapeBoxTypeTag { using InheritsFrom = std::tuple<BoxModel, PlumeShapeTypeTag>; };
} // end namespace TTag

//Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PlumeShapeTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PlumeShapeTypeTag> { using type = PlumeShapeProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PlumeShapeTypeTag>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = PlumeShapeSpatialParams<FVGridGeometry, Scalar>;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PlumeShapeTypeTag>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::BrineCO2<Scalar,
        Components::CO2<Scalar, GeneratedCO2Tables::CO2Tables>,
        Components::TabulatedComponent<Components::H2O<Scalar>>,
        FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/true>
    >;
};

// Use Moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::PlumeShapeTypeTag> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
