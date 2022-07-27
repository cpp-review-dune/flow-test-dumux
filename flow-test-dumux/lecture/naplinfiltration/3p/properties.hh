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
 * \brief Isothermal NAPL infiltration problem: L- or DNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
#ifndef DUMUX_NAPLINFILTRATIONPROPERTIES_3P_HH
#define DUMUX_NAPLINFILTRATIONPROPERTIES_3P_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p/model.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/3pimmiscible.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/air.hh>
#include "mymesitylene.hh"
#include <dumux/material/components/h2o.hh>

#include "../spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct InfiltrationThreePTypeTag { using InheritsFrom = std::tuple<ThreeP>; };
struct InfiltrationThreePBoxTypeTag { using InheritsFrom = std::tuple<InfiltrationThreePTypeTag, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::InfiltrationThreePTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::InfiltrationThreePTypeTag> { using type = InfiltrationThreePProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::InfiltrationThreePTypeTag>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Water = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using WettingFluid = FluidSystems::OnePLiquid<Scalar, Water>;
    using NonwettingFluid = FluidSystems::OnePLiquid<Scalar, Components::MyMesitylene<Scalar>>;
    using Gas = FluidSystems::OnePGas<Scalar, Components::Air<Scalar>>;
public:
    using type = FluidSystems::ThreePImmiscible<Scalar, WettingFluid, NonwettingFluid, Gas>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::InfiltrationThreePTypeTag>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InfiltrationSpatialParams<FVGridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
